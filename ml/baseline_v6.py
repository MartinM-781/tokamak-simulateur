"""Détecteurs baseline sur les features v6-diag : z-scores par diagnostic
(mono-canal), multi-canal, et régression logistique.

Usage :
    python ml/baseline_v6.py data/run_v6_diag02

Même protocole que la baseline v5 : split PAR TIR (70/30 stratifié, graine
fixe), statistiques robustes apprises sur les tirs sains du train, fenêtres
d'entraînement du classifieur = précurseur (HORIZON avant t_tq) contre
négatifs (tirs sains + fenêtres anciennes des tirs disruptifs), post-TQ exclu.
Écrit <run>/scores.csv (scores + z par feature pour l'analyse).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from features_v6 import FEATURES

# Détecteurs mono-diagnostic — z_saddle est l'équivalent du détecteur de mode
# verrouillé réel, z_mirnov celui de l'analyse de modes tournants.
CANAUX = {
    "z_mirnov": ["amp_m2", "amp_m3", "freq_m2_khz", "d_amp_m2"],
    "z_saddle": ["n1_amp", "n1_pente", "d_n1_amp"],
    "z_ece": ["te_coeur", "pente_te_coeur", "piquage", "te_ilot"],
    "z_bolo": ["bolo_mw", "pente_bolo"],
    "z_ip": ["ip_dev"],
    "z_nel": ["nel"],
    "z_multi": FEATURES,
}
HORIZON_S = 1.0   # fenêtre précurseur avant t_tq (échelles v6 : secondes)
MARGE_S = 2.5     # au-delà, les fenêtres disruptives servent de négatifs
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

    ref = df[(df["split"] == "train") & (df["label"] == "sain")]
    med = ref[FEATURES].median()
    mad = (ref[FEATURES] - med).abs().median() * 1.4826
    mad = mad.replace(0.0, 1e-12)
    z = ((df[FEATURES] - med) / mad).abs()
    for nom, cols in CANAUX.items():
        df[nom] = z[cols].max(axis=1)
    for col in FEATURES:
        df["z_" + col] = z[col]

    est_disruptif = df["label"] == "disruptif"
    avant_tq = df["t_s"] < df["t_tq_s"]
    precurseur = est_disruptif & avant_tq & (df["t_s"] >= df["t_tq_s"] - HORIZON_S)
    negatif = ~est_disruptif | (est_disruptif & (df["t_s"] < df["t_tq_s"] - MARGE_S))
    entrainable = (precurseur | negatif) & (df["split"] == "train")

    X_train = df.loc[entrainable, FEATURES].to_numpy()
    y_train = precurseur[entrainable].to_numpy()
    scaler = StandardScaler().fit(X_train)
    clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    clf.fit(scaler.transform(X_train), y_train)
    df["logistique"] = clf.predict_proba(scaler.transform(df[FEATURES].to_numpy()))[:, 1]

    cols = ["seed", "label", "split", "t_s", "phase", "t_tq_s",
            *CANAUX.keys(), "logistique", *("z_" + c for c in FEATURES)]
    dest = args.run_dir / "scores.csv"
    df[cols].to_csv(dest, index=False)

    print(f"{dest} : {len(df)} fenêtres scorées.")
    print(f"Split par tir : {len(train)} train / {len(test)} test "
          f"({int(y_train.sum())} fenêtres précurseur sur {len(y_train)} d'entraînement).")
    poids = pd.Series(clf.coef_[0], index=FEATURES).sort_values(key=abs, ascending=False)
    print("Poids logistiques (features standardisées) :")
    for nom, wgt in poids.items():
        print(f"  {nom:>15s} : {wgt:+.2f}")


if __name__ == "__main__":
    main()
