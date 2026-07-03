"""Features « physicien » sur les datasets diagnostics v6 (palier B).

Usage :
    python ml/features_v6.py data/run_v6_diag02 [--win-ms 10] [--step-ms 5]

Fait le travail d'une analyse de salle de contrôle sur les signaux BRUTS :
- décomposition en nombres de mode poloïdaux du réseau de Mirnov (DFT sur les
  8 bobines : amplitudes m=1/2/3, fréquence du mode m=2 par dérive de phase) ;
- amplitude n=1 des boucles à selle (détecteur de mode verrouillé) + pente ;
- ECE : cœur, pente du cœur, bord, piquage cœur/bord, canal traversé par
  l'îlot q=2 ;
- bolométrie (niveau + pente), interférométrie, écart d'Ip à sa médiane mobile.

Fenêtres glissantes horodatées FIN de fenêtre (aucun regard vers le futur).
Les colonnes de vérité terrain (phase, t_tq, label) ne servent qu'à
l'étiquetage et à l'évaluation en aval. Écrit <run>/features.csv.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view

FEATURES = [
    "amp_m1", "amp_m2", "amp_m3", "freq_m2_khz",
    "n1_amp", "n1_pente",
    "te_coeur", "pente_te_coeur", "te_bord", "piquage", "te_ilot",
    "bolo_mw", "pente_bolo", "nel", "ip_dev",
    "d_amp_m2", "d_n1_amp",
]


def charge_manifest(run_dir: Path) -> dict:
    with open(run_dir / "manifest.json", encoding="utf-8") as f:
        return json.load(f)


def fenetres(x: np.ndarray, win: int, hop: int) -> np.ndarray:
    return sliding_window_view(x, win)[::hop]


def pente(y_wins: np.ndarray, dt_s: float) -> np.ndarray:
    """Pente des moindres carrés (par seconde) de chaque fenêtre."""
    n = y_wins.shape[1]
    t = np.arange(n, dtype=float) * dt_s
    t_c = t - t.mean()
    return (y_wins @ t_c) / (t_c @ t_c)


def features_tir(df: pd.DataFrame, man: dict, win: int, hop: int) -> pd.DataFrame:
    dt_s = 1.0 / man["rate_hz"]
    cols = man["colonnes"]
    coils = df[[c for c in cols if c.startswith("mir_p")]].to_numpy(float)
    thetas = np.deg2rad(np.asarray(man["diagnostics"]["mirnov_thetas_deg"], float))
    sad = df[[c for c in cols if c.startswith("sad_t")]].to_numpy(float)
    ece = {c: df[c].to_numpy(float) for c in cols if c.startswith("ece_r")}
    bolo = df["bolo_MW"].to_numpy(float)
    nel = df["nel_1e19m2"].to_numpy(float)
    ip = df["ip_MA"].to_numpy(float)
    t = df["t_s"].to_numpy(float)
    phase = df["phase"].to_numpy(int)
    n = len(t)
    if n < win + hop:
        return pd.DataFrame()

    # --- Décomposition en modes : C_m(t) = moyenne_j coil_j·e^{-i m θ_j}
    expm = np.exp(-1j * np.outer(thetas, np.array([1.0, 2.0, 3.0])))  # (8, 3)
    C = coils @ expm / coils.shape[1]                                  # (n, 3)
    ampC = np.abs(C)
    # Fréquence du m=2 par dérive de phase entre échantillons consécutifs.
    dphi = np.angle(C[1:, 1] * np.conj(C[:-1, 1]))
    freq_hz = np.concatenate([[0.0], dphi]) / (2 * np.pi * dt_s)

    # --- Amplitude n=1 des boucles à selle (différences de boucles opposées)
    n1 = np.hypot(sad[:, 0] - sad[:, 2], sad[:, 1] - sad[:, 3]) / 2

    fin = np.arange(win - 1, n, hop)[: len(fenetres(t, win, hop))]
    w = lambda x: fenetres(x, win, hop)

    ip_med = np.median(w(ip), axis=1)
    # Baseline d'Ip : médiane mobile des médianes de fenêtres PRÉCÉDENTES (~1 s).
    n_hist = max(1, int(round(1.0 / (hop * dt_s))))
    ip_base = np.array(pd.Series(ip_med).shift(1).rolling(n_hist, min_periods=1).median())
    ip_base[0] = ip_med[0]

    out = pd.DataFrame(
        {
            "t_s": t[fin],
            "amp_m1": w(ampC[:, 0]).mean(axis=1),
            "amp_m2": w(ampC[:, 1]).mean(axis=1),
            "amp_m3": w(ampC[:, 2]).mean(axis=1),
            "freq_m2_khz": w(freq_hz).mean(axis=1) / 1e3,
            "n1_amp": w(n1).mean(axis=1),
            "n1_pente": pente(w(n1), dt_s),
            "te_coeur": w(ece["ece_r010"]).mean(axis=1),
            "pente_te_coeur": pente(w(ece["ece_r010"]), dt_s),
            "te_bord": w(ece["ece_r085"]).mean(axis=1),
            "te_ilot": w(ece["ece_r070"]).mean(axis=1),
            "bolo_mw": w(bolo).mean(axis=1),
            "pente_bolo": pente(w(bolo), dt_s),
            "nel": w(nel).mean(axis=1),
            "ip_dev": ip_med - ip_base,
            "phase": phase[fin],
        }
    )
    out["piquage"] = out["te_coeur"] / out["te_bord"].clip(lower=0.05)
    for col in ("amp_m2", "n1_amp"):
        out["d_" + col] = out[col].diff().fillna(0.0)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--win-ms", type=float, default=10.0)
    ap.add_argument("--step-ms", type=float, default=5.0)
    args = ap.parse_args()

    man = charge_manifest(args.run_dir)
    if man.get("schema") != "tokamak_v6_diag_run/1":
        raise SystemExit(f"schéma inattendu : {man.get('schema')} (attendu tokamak_v6_diag_run/1)")
    win = int(round(args.win_ms * 1e-3 * man["rate_hz"]))
    hop = int(round(args.step_ms * 1e-3 * man["rate_hz"]))

    morceaux = []
    for tir in man["tirs"]:
        df = pd.read_csv(args.run_dir / tir["fichier"], comment="#")
        feats = features_tir(df, man, win, hop)
        if feats.empty:
            continue
        feats.insert(0, "seed", tir["seed"])
        feats.insert(1, "label", tir["label"])
        feats["t_tq_s"] = tir["t_tq_s"] if tir["t_tq_s"] is not None else np.nan
        morceaux.append(feats)

    tout = pd.concat(morceaux, ignore_index=True)
    dest = args.run_dir / "features.csv"
    tout.to_csv(dest, index=False)
    print(f"{dest} : {len(tout)} fenêtres ({win} points à {man['rate_hz']} Hz, pas {hop}) "
          f"sur {tout['seed'].nunique()} tirs.")


if __name__ == "__main__":
    main()
