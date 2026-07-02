"""Extraction de features par fenêtres glissantes sur les tirs d'un run.

Usage :
    python ml/features.py data/run01 [--win-ms 20] [--step-ms 5]

Écrit <run>/features.csv : une ligne par fenêtre, horodatée à la FIN de la
fenêtre (aucun regard vers le futur). Features calculées sur les seuls canaux
MESURÉS (mirnov, te, ip) ; les colonnes de vérité terrain (label, phase, t_tq)
ne sont reportées que pour l'étiquetage et l'évaluation en aval.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view

SAMPLE_MS = 0.1  # période d'échantillonnage CSV (dt 0.05 ms x stride 2)


def charge_manifest(run_dir: Path) -> dict:
    with open(run_dir / "manifest.json", encoding="utf-8") as f:
        return json.load(f)


def fenetres(x: np.ndarray, win: int, hop: int) -> np.ndarray:
    """Vue (n_fenêtres, win) sans copie, une fenêtre tous les `hop` points."""
    return sliding_window_view(x, win)[::hop]


def pente(y_wins: np.ndarray, dt_ms: float) -> np.ndarray:
    """Pente des moindres carrés (par ms) de chaque fenêtre, vectorisée."""
    n = y_wins.shape[1]
    t = np.arange(n, dtype=float) * dt_ms
    t_c = t - t.mean()
    return (y_wins @ t_c) / (t_c @ t_c)


def features_tir(df: pd.DataFrame, win: int, hop: int) -> pd.DataFrame:
    mir = df["mirnov"].to_numpy(float)
    te = df["te"].to_numpy(float)
    ip = df["ip"].to_numpy(float)
    t = df["t_ms"].to_numpy(float)
    phase = df["phase"].to_numpy(int)
    if len(mir) < win + hop:
        return pd.DataFrame()

    w_mir = fenetres(mir, win, hop)
    w_te = fenetres(te, win, hop)
    w_ip = fenetres(ip, win, hop)
    fin = np.arange(win - 1, len(mir), hop)[: len(w_mir)]  # index de fin de fenêtre

    rms_mirnov = np.sqrt((w_mir**2).mean(axis=1))
    # Fréquence par passages à zéro : n_croisements / (2 * durée) — en kHz car 1/ms = kHz.
    croisements = (np.diff(np.signbit(w_mir), axis=1) != 0).sum(axis=1)
    freq_khz = croisements / (2.0 * win * SAMPLE_MS)
    te_slope = pente(w_te, SAMPLE_MS)
    # Écart d'ip à sa médiane mobile : médiane de la fenêtre vs médiane des
    # 100 ms qui PRÉCÈDENT la fenêtre (baseline non contaminée par la fenêtre).
    ip_med = pd.Series(ip).rolling(1000, min_periods=win).median().to_numpy()
    debut = fin - (win - 1)
    baseline = ip_med[np.maximum(debut - 1, win - 1)]
    ip_dev = np.median(w_ip, axis=1) - baseline

    out = pd.DataFrame(
        {
            "t_ms": t[fin],
            "rms_mirnov": rms_mirnov,
            "freq_khz": freq_khz,
            "te_slope": te_slope,
            "ip_dev": ip_dev,
            "phase": phase[fin],
        }
    )
    # Dérivées inter-fenêtres (par pas de hop) — première fenêtre : 0.
    for col in ("rms_mirnov", "freq_khz", "te_slope", "ip_dev"):
        out["d_" + col] = out[col].diff().fillna(0.0)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--win-ms", type=float, default=20.0)
    ap.add_argument("--step-ms", type=float, default=5.0)
    args = ap.parse_args()

    win = int(round(args.win_ms / SAMPLE_MS))
    hop = int(round(args.step_ms / SAMPLE_MS))
    manifest = charge_manifest(args.run_dir)

    morceaux = []
    for tir in manifest["tirs"]:
        df = pd.read_csv(args.run_dir / tir["fichier"], comment="#")
        feats = features_tir(df, win, hop)
        if feats.empty:
            continue
        feats.insert(0, "seed", tir["seed"])
        feats.insert(1, "label", tir["label"])
        feats["t_tq_ms"] = tir["t_tq_ms"] if tir["t_tq_ms"] is not None else np.nan
        morceaux.append(feats)

    tout = pd.concat(morceaux, ignore_index=True)
    dest = args.run_dir / "features.csv"
    tout.to_csv(dest, index=False)
    print(
        f"{dest} : {len(tout)} fenêtres ({win} points, pas {hop}) "
        f"sur {tout['seed'].nunique()} tirs."
    )


if __name__ == "__main__":
    main()
