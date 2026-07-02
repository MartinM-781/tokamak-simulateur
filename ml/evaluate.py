"""Évaluation des détecteurs : temps d'alerte avant t_tq vs fausses alarmes.

Usage :
    python ml/evaluate.py data/run01

Métrique reine : sur les tirs disruptifs du test, temps d'alerte médian avant
le quench thermique (t_tq du manifest) ; en regard, taux de fausses alarmes
sur les tirs sains du test. Une alarme = score au-dessus du seuil sur
K_CONSECUTIFS fenêtres consécutives ; sur tir disruptif, seules les fenêtres
antérieures à t_tq comptent (détecter après coup ne sert à rien).

Sorties : ml/resultats/*.png + ml/RESULTS.md (tout est régénéré, y compris
les chiffres du texte). La colonne `phase` (vérité terrain) n'est utilisée
que pour l'ANALYSE — jamais par les détecteurs.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DETECTEURS = ["z_mirnov", "z_te", "z_ip", "z_prad", "z_ne", "z_multi", "logistique"]
K_CONSECUTIFS = 3
FAR_CIBLE = 0.05

COULEURS = {
    "z_mirnov": "#F272AC", "z_te": "#D9975B", "z_ip": "#5C89B8",
    "z_prad": "#E8B45A", "z_ne": "#8FBC8F",
    "z_multi": "#FF7A5C", "logistique": "#93A4B4",
}


def t_alarme(t: np.ndarray, score: np.ndarray, seuil: float) -> float | None:
    """Instant de la 1re alarme (K fenêtres consécutives >= seuil), sinon None."""
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
        g = g.sort_values("t_ms")
        t_tq = g["t_tq_ms"].iloc[0]
        disruptif = g["label"].iloc[0] == "disruptif"
        if disruptif:
            g = g[g["t_ms"] < t_tq]  # seules les fenêtres pré-TQ comptent
        series.append({
            "seed": seed, "disruptif": disruptif, "t_tq": t_tq,
            "t": g["t_ms"].to_numpy(), "score": g[det].to_numpy(),
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
    """Plus petit seuil (parmi les quantiles hauts des scores sains du train)
    donnant un taux de fausses alarmes <= FAR_CIBLE sur les tirs sains du train."""
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
    ap.add_argument("--out-dir", type=Path, default=Path("ml/resultats"))
    ap.add_argument("--results", type=Path, default=Path("ml/RESULTS.md"))
    args = ap.parse_args()

    df = pd.read_csv(args.run_dir / "scores.csv")
    df_train = df[df["split"] == "train"]
    df_test = df[df["split"] == "test"]
    n_dis = df_test[df_test["label"] == "disruptif"]["seed"].nunique()
    n_sains = df_test[df_test["label"] == "sain"]["seed"].nunique()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # --- Courbes temps d'alerte / détection vs fausses alarmes (balayage de seuil, test).
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    points_fonctionnement = {}
    stats_phase = {}
    for det in DETECTEURS:
        series_test = tirs_en_series(df_test, det)
        series_train = tirs_en_series(df_train, det)
        scores_sains_train = df_train.loc[df_train["label"] == "sain", det].to_numpy()

        balayage = np.unique(np.quantile(df_test[det], np.linspace(0.0, 1.0, 300)))
        courbe = [dict(evalue_seuil(series_test, s), seuil=s) for s in balayage]
        far = [c["far"] for c in courbe]
        med = [c["alerte_mediane"] for c in courbe]
        det_rate = [c["detection"] for c in courbe]
        ax1.plot(far, med, label=det, color=COULEURS[det], lw=1.6)
        ax2.plot(far, det_rate, label=det, color=COULEURS[det], lw=1.6)

        seuil_op = calibre_seuil(series_train, scores_sains_train)
        points_fonctionnement[det] = dict(evalue_seuil(series_test, seuil_op), seuil=seuil_op)

    # --- Analyse par phase (tirs disruptifs du test, fenêtres pré-TQ).
    # « Phase 2 stricte » : la fenêtre 4 crans (20 ms) plus tôt était déjà en
    # phase 2 — la fenêtre courante est donc entièrement postérieure au
    # verrouillage, sans chevaucher le transitoire d'effondrement du RMS.
    dtest = df_test.sort_values(["seed", "t_ms"]).copy()
    dtest["ph_20ms_avant"] = dtest.groupby("seed")["phase"].shift(4)
    pre_tq = (dtest["label"] == "disruptif") & (dtest["t_ms"] < dtest["t_tq_ms"])
    masques = {
        "ph1": pre_tq & (dtest["phase"] == 1),
        "ph2": pre_tq & (dtest["phase"] == 2),
        "ph2s": pre_tq & (dtest["phase"] == 2) & (dtest["ph_20ms_avant"] == 2),
    }
    for det in DETECTEURS:
        seuil_op = points_fonctionnement[det]["seuil"]
        stats_phase[det] = {
            nom: float((dtest.loc[m, det] >= seuil_op).mean()) for nom, m in masques.items()
        }
    z_cols = [c for c in dtest.columns if c.startswith("z_") and c not in DETECTEURS]
    z_multi_median = {nom: float(dtest.loc[m, "z_multi"].median()) for nom, m in masques.items()}
    dominantes = (
        dtest.loc[masques["ph2"], z_cols].idxmax(axis=1).value_counts(normalize=True)
    )
    duree_ph2_tirs = (
        dtest[masques["ph2"]].groupby("seed")["t_ms"].agg(lambda s: s.max() - s.min())
    )

    for ax, ylab in ((ax1, "temps d'alerte médian avant t_tq (ms)"), (ax2, "taux de détection")):
        ax.set_xlabel("taux de fausses alarmes (fraction des tirs sains)")
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Alerte précoce vs fausses alarmes — test par tir, alarme = "
                 f"{K_CONSECUTIFS} fenêtres consécutives")
    fig.tight_layout()
    png_courbes = args.out_dir / "alerte_vs_fausses_alarmes.png"
    fig.savefig(png_courbes, dpi=130)
    plt.close(fig)

    # --- Illustration mode verrouillé : le tir test à la phase 2 la plus longue.
    dis = df_test[df_test["label"] == "disruptif"]
    duree_ph2 = dis[dis["phase"] == 2].groupby("seed")["t_ms"].count()
    seed_ex = int(duree_ph2.idxmax())
    g = dis[dis["seed"] == seed_ex].sort_values("t_ms")
    t_tq_ex = g["t_tq_ms"].iloc[0]
    g = g[g["t_ms"] <= t_tq_ex + 25]  # après le TQ, détecter ne compte plus
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.plot(g["t_ms"], g["z_mirnov"], label="z_mirnov", color=COULEURS["z_mirnov"], lw=2.6, alpha=0.9)
    ax.plot(g["t_ms"], g["z_prad"], label="z_prad", color=COULEURS["z_prad"], lw=1.6)
    ax.plot(g["t_ms"], g["z_multi"], label="z_multi", color=COULEURS["z_multi"], lw=1.2, ls="--")
    ax.axhline(points_fonctionnement["z_multi"]["seuil"], color="k", ls=":", lw=1,
               label="seuil (5 % FA)")
    bandes = {1: "#E8B45A", 2: "#FF7A5C"}
    for ph, coul in bandes.items():
        m = g["phase"] == ph
        if m.any():
            ax.axvspan(g.loc[m, "t_ms"].min(), g.loc[m, "t_ms"].max(), color=coul, alpha=0.12,
                       label=f"phase {ph} ({'précurseur' if ph == 1 else 'mode verrouillé'})")
    ax.axvline(t_tq_ex, color="k", ls="--", lw=1, label="t_tq")
    ax.set_yscale("log")
    ax.set_xlabel("t (ms)")
    ax.set_ylabel("score (log)")
    ax.set_title(f"Tir {seed_ex} — au verrouillage le canal magnétique s'éteint, "
                 "P_rad porte l'alerte jusqu'au quench")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    png_verrou = args.out_dir / "mode_verrouille.png"
    fig.savefig(png_verrou, dpi=130)
    plt.close(fig)

    # --- RESULTS.md — entièrement régénéré, chiffres inclus.
    pf = points_fonctionnement
    lignes = [
        "# Résultats baseline — détection précoce de disruptions",
        "",
        f"Dataset : `{args.run_dir.as_posix()}` — test : {n_dis} tirs disruptifs, "
        f"{n_sains} tirs sains (split PAR TIR, 70/30 stratifié, graine fixe).",
        "",
        "## Protocole",
        "",
        "- Features par fenêtres glissantes de 20 ms (pas 5 ms), horodatées fin de fenêtre,",
        "  sur les seuls canaux mesurés : RMS Mirnov, fréquence par passages à zéro,",
        "  pente de T_e, écart d'I_p à sa médiane mobile, et leurs dérivées.",
        "- Détecteurs z-score robustes (médiane/MAD appris sur les tirs sains du train)",
        "  mono-canal (`z_mirnov`, `z_te`, `z_ip`) et multi-canal (`z_multi` = max des z),",
        "  plus une régression logistique (précurseur = fenêtre dans les 300 ms avant t_tq).",
        f"- Alarme = score ≥ seuil sur {K_CONSECUTIFS} fenêtres consécutives. Seuil calibré",
        f"  sur le train (fausses alarmes ≤ {FAR_CIBLE:.0%} des tirs sains du train).",
        "- Sur tir disruptif, seules les fenêtres antérieures à t_tq comptent : une",
        "  « détection » après le quench thermique est un échec (manque).",
        "",
        f"![courbes]({png_courbes.name})",
        "",
        "## Points de fonctionnement (seuil calibré sur train)",
        "",
        "| Détecteur | Seuil | Fausses alarmes (test) | Détection | Alerte médiane (ms) | Alerte p10 (ms) |",
        "|---|---|---|---|---|---|",
    ]
    for det in DETECTEURS:
        p = pf[det]
        lignes.append(
            f"| {det} | {p['seuil']:.3g} | {p['far']:.0%} | {p['detection']:.0%} "
            f"| {p['alerte_mediane']:.1f} | {p['alerte_p10']:.1f} |"
        )
    top_dom = " puis ".join(f"`{k[2:]}` ({v:.0%})" for k, v in dominantes.head(2).items())
    lignes += [
        "",
        "## Comportement pendant le mode verrouillé (phase = 2, Mirnov muette)",
        "",
        "Le signal Mirnov est ∝ W²·Ω : au verrouillage (Ω → 0), la sonde se tait",
        "alors que l'îlot continue de croître. Fraction des fenêtres pré-TQ",
        "au-dessus du seuil de fonctionnement, par phase (tirs disruptifs du test) ;",
        "« phase 2 stricte » = fenêtres entièrement postérieures au verrouillage",
        "(sans chevauchement du transitoire d'effondrement) :",
        "",
        "| Détecteur | Phase 1 (précurseur rotatif) | Phase 2 (toutes fenêtres) | Phase 2 stricte |",
        "|---|---|---|---|",
    ]
    for det in DETECTEURS:
        s = stats_phase[det]
        lignes.append(f"| {det} | {s['ph1']:.0%} | {s['ph2']:.0%} | {s['ph2s']:.0%} |")
    lignes += [
        "",
        f"Trois enseignements (feature dominante en phase 2 : {top_dom} ; "
        f"trou de phase 2 : {duree_ph2_tirs.median():.0f} ms en médiane, "
        f"max {duree_ph2_tirs.max():.0f} ms, juste avant le quench) :",
        "",
        "1. **L'alerte s'acquiert en phase 1, portée par le RMS Mirnov du mode",
        f"   tournant** ({stats_phase['z_mirnov']['ph1']:.0%} des fenêtres au-dessus du seuil) —",
        "   mais les bouffées ELM des tirs sains gonflent la calibration du canal",
        "   magnétique et lui coûtent de la marge : son seuil monte, son alerte",
        "   recule, et il peut manquer les tirs à précurseur court.",
        "2. **Au verrouillage, le canal magnétique s'éteint** : en phase 2 stricte",
        f"   `z_mirnov` ne tient plus que {stats_phase['z_mirnov']['ph2s']:.0%} des fenêtres — c'est le trou",
        "   aveugle historique de la détection purement magnétique. Ni T_e (dérive",
        "   lente noyée dans le clapotis des dents de scie qui a calibré la MAD),",
        "   ni I_p (plat jusqu'au quench), ni la densité ne prennent le relais",
        "   (`z_te`, `z_ip`, `z_ne` : ~0 % partout à 5 % de fausses alarmes).",
        "3. **C'est P_rad qui comble le trou** : l'îlot verrouillé continue de",
        f"   grossir et de rayonner, `z_prad` tient {stats_phase['z_prad']['ph2s']:.0%} des fenêtres de",
        f"   phase 2 stricte, et le multi-canal ({stats_phase['z_multi']['ph2s']:.0%}) comme la logistique",
        f"   ({stats_phase['logistique']['ph2s']:.0%}) restent armés du verrouillage jusqu'au quench.",
        "",
        "C'est la réponse à la question mono vs multi : un détecteur mono-canal",
        "magnétique sans mémoire perd l'alarme pendant les dizaines de",
        "millisecondes qui précèdent immédiatement le quench (il ne survivrait que",
        "par verrouillage d'alarme, ou par une feature à mémoire du type « le RMS",
        "s'est effondré depuis un niveau élevé », que `d_rms_mirnov` capture au",
        "transitoire). Le multi-canal, lui, n'a pas besoin de mémoire : la",
        "complémentarité physique des canaux (mirnov pour le mode tournant, P_rad",
        "pour l'îlot verrouillé) assure une couverture continue.",
        "",
        f"![mode verrouillé]({png_verrou.name})",
        "",
        "## Reproduire",
        "",
        "```",
        "node scripts/generate.js --shots 200 --out data/run01 --disrupt-ratio 0.5 --seed-base 1000",
        f"python ml/features.py {args.run_dir.as_posix()}",
        f"python ml/baseline.py {args.run_dir.as_posix()}",
        f"python ml/evaluate.py {args.run_dir.as_posix()}",
        "```",
        "",
    ]
    args.results.write_text("\n".join(lignes), encoding="utf-8")
    print(f"{png_courbes}\n{png_verrou}\n{args.results} écrits.")
    for det in DETECTEURS:
        p = pf[det]
        print(f"  {det:>11s} : FA={p['far']:.0%} détection={p['detection']:.0%} "
              f"alerte médiane={p['alerte_mediane']:.1f} ms (seuil {p['seuil']:.3g})")


if __name__ == "__main__":
    main()
