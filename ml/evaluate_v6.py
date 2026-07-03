"""Évaluation des détecteurs v6-diag : temps d'alerte avant t_tq vs fausses
alarmes, et comportement pendant le mode verrouillé (Mirnov vs boucles à selle).

Usage :
    python ml/evaluate_v6.py data/run_v6_diag02

Alarme = score ≥ seuil sur K_CONSECUTIFS fenêtres consécutives ; seuil calibré
sur les tirs sains du train (fausses alarmes ≤ FAR_CIBLE) ; sur tir disruptif,
seules les fenêtres antérieures à t_tq comptent. Sorties : ml/resultats_v6/*.png
+ ml/RESULTS_V6.md, entièrement régénérés (chiffres compris). La colonne
`phase` (vérité terrain) ne sert qu'à l'ANALYSE, jamais aux détecteurs.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DETECTEURS = ["z_mirnov", "z_saddle", "z_ece", "z_bolo", "z_ip", "z_nel",
              "z_multi", "logistique"]
K_CONSECUTIFS = 3
FAR_CIBLE = 0.05

COULEURS = {
    "z_mirnov": "#F272AC", "z_saddle": "#5C89B8", "z_ece": "#D9975B",
    "z_bolo": "#E8B45A", "z_ip": "#8FBC8F", "z_nel": "#B39DDB",
    "z_multi": "#FF7A5C", "logistique": "#93A4B4",
}


def t_alarme(t: np.ndarray, score: np.ndarray, seuil: float) -> float | None:
    au_dessus = score >= seuil
    consec = 0
    for i in range(len(au_dessus)):
        consec = consec + 1 if au_dessus[i] else 0
        if consec >= K_CONSECUTIFS:
            return float(t[i])
    return None


def tirs_en_series(df: pd.DataFrame, det: str) -> list[dict]:
    series = []
    for seed, g in df.groupby("seed"):
        g = g.sort_values("t_s")
        t_tq = g["t_tq_s"].iloc[0]
        disruptif = g["label"].iloc[0] == "disruptif"
        if disruptif:
            g = g[g["t_s"] < t_tq]
        series.append({
            "seed": seed, "disruptif": disruptif, "t_tq": t_tq,
            "t": g["t_s"].to_numpy(), "score": g[det].to_numpy(),
        })
    return series


def evalue_seuil(series: list[dict], seuil: float) -> dict:
    alertes, manques, fausses, n_sains = [], 0, 0, 0
    for s in series:
        ta = t_alarme(s["t"], s["score"], seuil)
        if s["disruptif"]:
            if ta is None:
                manques += 1
            else:
                alertes.append(s["t_tq"] - ta)
        else:
            n_sains += 1
            if ta is not None:
                fausses += 1
    n_dis = manques + len(alertes)
    return {
        "far": fausses / max(n_sains, 1),
        "detection": len(alertes) / max(n_dis, 1),
        "alerte_mediane": float(np.median(alertes)) if alertes else np.nan,
        "alerte_p10": float(np.percentile(alertes, 10)) if alertes else np.nan,
    }


def calibre_seuil(series_train: list[dict], scores_train_sains: np.ndarray) -> float:
    candidats = np.unique(np.quantile(scores_train_sains, np.linspace(0.90, 1.0, 200)))
    sains = [s for s in series_train if not s["disruptif"]]
    for seuil in candidats:
        far = sum(t_alarme(s["t"], s["score"], seuil) is not None for s in sains) / max(len(sains), 1)
        if far <= FAR_CIBLE:
            return float(seuil)
    return float(candidats[-1])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out-dir", type=Path, default=Path("ml/resultats_v6"))
    ap.add_argument("--results", type=Path, default=Path("ml/RESULTS_V6.md"))
    args = ap.parse_args()

    df = pd.read_csv(args.run_dir / "scores.csv")
    df_train = df[df["split"] == "train"]
    df_test = df[df["split"] == "test"]
    n_dis = df_test[df_test["label"] == "disruptif"]["seed"].nunique()
    n_sains = df_test[df_test["label"] == "sain"]["seed"].nunique()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    pf = {}
    for det in DETECTEURS:
        series_test = tirs_en_series(df_test, det)
        series_train = tirs_en_series(df_train, det)
        scores_sains_train = df_train.loc[df_train["label"] == "sain", det].to_numpy()

        balayage = np.unique(np.quantile(df_test[det], np.linspace(0.0, 1.0, 300)))
        courbe = [dict(evalue_seuil(series_test, s), seuil=s) for s in balayage]
        ax1.plot([c["far"] for c in courbe], [c["alerte_mediane"] for c in courbe],
                 label=det, color=COULEURS[det], lw=1.6)
        ax2.plot([c["far"] for c in courbe], [c["detection"] for c in courbe],
                 label=det, color=COULEURS[det], lw=1.6)

        seuil_op = calibre_seuil(series_train, scores_sains_train)
        pf[det] = dict(evalue_seuil(series_test, seuil_op), seuil=seuil_op)

    for ax, ylab in ((ax1, "temps d'alerte médian avant t_tq (s)"), (ax2, "taux de détection")):
        ax.set_xlabel("taux de fausses alarmes (fraction des tirs sains)")
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle(f"v6-diag — alerte précoce vs fausses alarmes (alarme = {K_CONSECUTIFS} fenêtres consécutives)")
    fig.tight_layout()
    png_courbes = args.out_dir / "alerte_vs_fausses_alarmes.png"
    fig.savefig(png_courbes, dpi=130)
    plt.close(fig)

    # --- Analyse par phase (fenêtres pré-TQ des tirs disruptifs du test)
    dtest = df_test.sort_values(["seed", "t_s"]).copy()
    dtest["ph_avant"] = dtest.groupby("seed")["phase"].shift(4)
    pre_tq = (dtest["label"] == "disruptif") & (dtest["t_s"] < dtest["t_tq_s"])
    masques = {
        "ph1": pre_tq & (dtest["phase"] == 1),
        "ph2": pre_tq & (dtest["phase"] == 2),
        "ph2s": pre_tq & (dtest["phase"] == 2) & (dtest["ph_avant"] == 2),
    }
    stats_phase = {
        det: {nom: float((dtest.loc[m, det] >= pf[det]["seuil"]).mean()) for nom, m in masques.items()}
        for det in DETECTEURS
    }
    z_cols = [c for c in dtest.columns if c.startswith("z_") and c not in DETECTEURS]
    dominantes = dtest.loc[masques["ph2s"], z_cols].idxmax(axis=1).value_counts(normalize=True)
    duree_ph2 = dtest[masques["ph2"]].groupby("seed")["t_s"].agg(lambda s: s.max() - s.min())

    # --- Illustration : le tir test à la phase verrouillée la plus longue
    dis = dtest[dtest["label"] == "disruptif"]
    seed_ex = int(dis[dis["phase"] == 2].groupby("seed")["t_s"].count().idxmax())
    g = dis[dis["seed"] == seed_ex].sort_values("t_s")
    t_tq_ex = g["t_tq_s"].iloc[0]
    g = g[g["t_s"] <= t_tq_ex + 0.2]
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.plot(g["t_s"], g["z_mirnov"], label="z_mirnov (modes tournants)", color=COULEURS["z_mirnov"], lw=2.4, alpha=0.9)
    ax.plot(g["t_s"], g["z_saddle"], label="z_saddle (mode verrouillé)", color=COULEURS["z_saddle"], lw=1.8)
    ax.plot(g["t_s"], g["z_multi"], label="z_multi", color=COULEURS["z_multi"], lw=1.1, ls="--")
    ax.axhline(pf["z_multi"]["seuil"], color="k", ls=":", lw=1, label="seuil (5 % FA)")
    for ph, coul, nom in ((1, "#E8B45A", "précurseur"), (2, "#FF7A5C", "mode verrouillé")):
        m = g["phase"] == ph
        if m.any():
            ax.axvspan(g.loc[m, "t_s"].min(), g.loc[m, "t_s"].max(), color=coul, alpha=0.12,
                       label=f"phase {ph} ({nom})")
    ax.axvline(t_tq_ex, color="k", ls="--", lw=1, label="t_tq")
    ax.set_yscale("log")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("score (log)")
    ax.set_title(f"Tir {seed_ex} — au verrouillage, la selle prend le relais des Mirnov")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    png_verrou = args.out_dir / "mode_verrouille.png"
    fig.savefig(png_verrou, dpi=130)
    plt.close(fig)

    # --- RESULTS_V6.md
    top_dom = " puis ".join(f"`{k[2:]}` ({v:.0%})" for k, v in dominantes.head(2).items())
    lignes = [
        "# Résultats baseline v6-diag — détection précoce sur signaux de salle de contrôle",
        "",
        f"Dataset : `{args.run_dir.as_posix()}` — test : {n_dis} tirs disruptifs, "
        f"{n_sains} tirs sains (split PAR TIR, 70/30 stratifié, graine fixe).",
        "",
        "## Protocole",
        "",
        "- Features par fenêtres de 10 ms (pas 5 ms) sur les 22 canaux bruts :",
        "  décomposition en modes m=1/2/3 du réseau de Mirnov (DFT sur les 8",
        "  bobines) + fréquence du m=2 par dérive de phase, amplitude n=1 des",
        "  boucles à selle (+ pente), ECE (cœur, pente, bord, piquage, canal q=2),",
        "  bolométrie, interférométrie, écart d'Ip, dérivées inter-fenêtres.",
        "- Détecteurs z-score robustes (médiane/MAD sur tirs sains du train)",
        "  par diagnostic (`z_mirnov`, `z_saddle`, `z_ece`, `z_bolo`, `z_ip`,",
        "  `z_nel`), multi-diagnostic (`z_multi`), et régression logistique",
        f"  (précurseur = {1.0:.0f} s avant t_tq).",
        f"- Alarme = score ≥ seuil sur {K_CONSECUTIFS} fenêtres consécutives ;",
        f"  seuil calibré sur le train (FA ≤ {FAR_CIBLE:.0%} des tirs sains).",
        "- Sur tir disruptif, seules les fenêtres antérieures à t_tq comptent.",
        "",
        f"![courbes]({png_courbes.name})",
        "",
        "## Points de fonctionnement (seuil calibré sur train)",
        "",
        "| Détecteur | Seuil | Fausses alarmes (test) | Détection | Alerte médiane (s) | Alerte p10 (s) |",
        "|---|---|---|---|---|---|",
    ]
    for det in DETECTEURS:
        p = pf[det]
        lignes.append(
            f"| {det} | {p['seuil']:.10g} | {p['far']:.0%} | {p['detection']:.0%} "
            f"| {p['alerte_mediane']:.2f} | {p['alerte_p10']:.2f} |"
        )
    depasse = [d for d in DETECTEURS if pf[d]["far"] > FAR_CIBLE]
    if depasse:
        lignes += ["", "Nota : " + ", ".join(f"`{d}` ({pf[d]['far']:.0%})" for d in depasse) +
                   f" dépasse{'nt' if len(depasse) > 1 else ''} la cible de {FAR_CIBLE:.0%} hors échantillon"
                   f" ({n_sains} tirs sains au test : incertitude binomiale large)."]
    lignes += [
        "",
        "## Mode verrouillé : Mirnov vs boucles à selle (phase = 2)",
        "",
        "Fraction des fenêtres pré-TQ au-dessus du seuil, par phase ;",
        "« phase 2 stricte » = fenêtres entièrement postérieures au verrouillage :",
        "",
        "| Détecteur | Phase 1 (mode tournant) | Phase 2 (toutes) | Phase 2 stricte |",
        "|---|---|---|---|",
    ]
    for det in DETECTEURS:
        s = stats_phase[det]
        lignes.append(f"| {det} | {s['ph1']:.0%} | {s['ph2']:.0%} | {s['ph2s']:.0%} |")
    lignes += [
        "",
        f"Feature dominante en phase 2 stricte : {top_dom}. Durée de la phase",
        f"verrouillée : médiane {duree_ph2.median():.2f} s, max {duree_ph2.max():.2f} s.",
        "",
        f"Lecture : en phase verrouillée stricte, `z_mirnov` tient "
        f"{stats_phase['z_mirnov']['ph2s']:.0%} des fenêtres (le dB/dt s'éteint avec la rotation)"
        f" tandis que `z_saddle` en tient {stats_phase['z_saddle']['ph2s']:.0%} : le champ radial",
        "statique de l'îlot verrouillé continue de croître, exactement ce que",
        "mesure un détecteur de mode verrouillé réel. C'est la réponse",
        "instrumentale au trou aveugle magnétique : il ne se comble pas par de",
        "la mémoire d'alarme mais par le BON capteur. Le multi-diagnostic",
        f"(`z_multi` : {stats_phase['z_multi']['ph2s']:.0%}) et la logistique"
        f" ({stats_phase['logistique']['ph2s']:.0%}) en héritent.",
        "",
        f"![mode verrouillé]({png_verrou.name})",
        "",
        "## Reproduire",
        "",
        "```",
        "node scripts/generate_v6_diags.js --shots 100 --out data/run_v6_diag02 --disrupt-ratio 0.5 --seed-base 3000",
        f"python ml/features_v6.py {args.run_dir.as_posix()}",
        f"python ml/baseline_v6.py {args.run_dir.as_posix()}",
        f"python ml/evaluate_v6.py {args.run_dir.as_posix()}",
        "```",
        "",
    ]
    args.results.write_text("\n".join(lignes), encoding="utf-8")
    print(f"{png_courbes}\n{png_verrou}\n{args.results} écrits.")
    for det in DETECTEURS:
        p = pf[det]
        print(f"  {det:>11s} : FA={p['far']:.0%} détection={p['detection']:.0%} "
              f"alerte médiane={p['alerte_mediane']:.2f} s (seuil {p['seuil']:.3g})")


if __name__ == "__main__":
    main()
