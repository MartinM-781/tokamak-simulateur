"""Détecteurs baseline : z-score robuste mono/multi-canal + régression logistique.

Usage :
    python ml/baseline.py data/run01

Lit <run>/features.csv, écrit <run>/scores.csv (un score par fenêtre et par
détecteur). Split train/test PAR TIR (jamais par fenêtre), stratifié par
label, graine fixe — les fenêtres d'un même tir ne sont jamais réparties
entre train et test.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

FEATURES = [
    "rms_mirnov", "freq_khz", "te_slope", "ip_dev",
    "prad_mean", "prad_slope", "ne_mean", "ne_slope",
    "d_rms_mirnov", "d_freq_khz", "d_te_slope", "d_ip_dev",
]
# Jeux de features par détecteur z-score : mono-canal vs multi-canal.
CANAUX = {
    "z_mirnov": ["rms_mirnov", "freq_khz", "d_rms_mirnov", "d_freq_khz"],
    "z_te": ["te_slope", "d_te_slope"],
    "z_ip": ["ip_dev", "d_ip_dev"],
    "z_prad": ["prad_mean", "prad_slope"],
    "z_ne": ["ne_mean", "ne_slope"],
    "z_multi": FEATURES,
}
# Étiquetage des fenêtres pour le classifieur : précurseur = fenêtre finissant
# dans les HORIZON_MS avant le TQ ; les fenêtres d'un tir disruptif plus
# anciennes que MARGE_MS avant le TQ servent de négatifs ; entre les deux,
# zone grise exclue de l'entraînement ; tout ce qui suit le TQ est exclu
# (détecter une disruption déjà en cours n'a pas d'intérêt).
HORIZON_MS = 300.0
MARGE_MS = 600.0
SEED_SPLIT = 0
PART_TEST = 0.3


def split_par_tir(df: pd.DataFrame) -> tuple[set[int], set[int]]:
    rng = np.random.default_rng(SEED_SPLIT)
    train: set[int] = set()
    test: set[int] = set()
    for label in ("disruptif", "sain"):
        seeds = sorted(df.loc[df["label"] == label, "seed"].unique())
        rng.shuffle(seeds)
        n_test = int(round(len(seeds) * PART_TEST))
        test.update(seeds[:n_test])
        train.update(seeds[n_test:])
    return train, test


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    args = ap.parse_args()

    df = pd.read_csv(args.run_dir / "features.csv")
    train, test = split_par_tir(df)
    df["split"] = np.where(df["seed"].isin(list(test)), "test", "train")

    # --- Détecteurs z-score : stats robustes apprises sur les tirs SAINS du train.
    ref = df[(df["split"] == "train") & (df["label"] == "sain")]
    med = ref[FEATURES].median()
    mad = (ref[FEATURES] - med).abs().median() * 1.4826
    mad = mad.replace(0.0, 1e-12)
    z = ((df[FEATURES] - med) / mad).abs()
    for nom, cols in CANAUX.items():
        df[nom] = z[cols].max(axis=1)
    # z par feature, exporté pour l'analyse fine (evaluate.py : quelle feature
    # domine pendant le mode verrouillé ?).
    for col in FEATURES:
        df["z_" + col] = z[col]

    # --- Classifieur logistique sur fenêtres étiquetées, split par tir.
    est_disruptif = df["label"] == "disruptif"
    avant_tq = df["t_ms"] < df["t_tq_ms"]
    precurseur = est_disruptif & avant_tq & (df["t_ms"] >= df["t_tq_ms"] - HORIZON_MS)
    negatif = ~est_disruptif | (est_disruptif & (df["t_ms"] < df["t_tq_ms"] - MARGE_MS))
    entrainable = (precurseur | negatif) & (df["split"] == "train")

    X_train = df.loc[entrainable, FEATURES].to_numpy()
    y_train = precurseur[entrainable].to_numpy()
    scaler = StandardScaler().fit(X_train)
    clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    clf.fit(scaler.transform(X_train), y_train)
    df["logistique"] = clf.predict_proba(scaler.transform(df[FEATURES].to_numpy()))[:, 1]

    cols = ["seed", "label", "split", "t_ms", "phase", "t_tq_ms",
            *CANAUX.keys(), "logistique", *("z_" + c for c in FEATURES)]
    dest = args.run_dir / "scores.csv"
    df[cols].to_csv(dest, index=False)

    n_pos = int(y_train.sum())
    print(f"{dest} : {len(df)} fenêtres scorées.")
    print(f"Split par tir : {len(train)} train / {len(test)} test "
          f"(fenêtres d'entraînement du classifieur : {len(y_train)}, dont {n_pos} précurseurs).")
    poids = pd.Series(clf.coef_[0], index=FEATURES).sort_values(key=abs, ascending=False)
    print("Poids logistiques (features standardisées) :")
    for nom, w in poids.items():
        print(f"  {nom:>14s} : {w:+.2f}")


if __name__ == "__main__":
    main()
