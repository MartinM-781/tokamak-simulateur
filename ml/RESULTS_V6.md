# Résultats baseline v6-diag — détection précoce sur signaux de salle de contrôle

Dataset : `data/run_v6_diag02` — test : 15 tirs disruptifs, 15 tirs sains (split PAR TIR, 70/30 stratifié, graine fixe).

## Protocole

- Features par fenêtres de 10 ms (pas 5 ms) sur les 22 canaux bruts :
  décomposition en modes m=1/2/3 du réseau de Mirnov (DFT sur les 8
  bobines) + fréquence du m=2 par dérive de phase, amplitude n=1 des
  boucles à selle (+ pente), ECE (cœur, pente, bord, piquage, canal q=2),
  bolométrie, interférométrie, écart d'Ip, dérivées inter-fenêtres.
- Détecteurs z-score robustes (médiane/MAD sur tirs sains du train)
  par diagnostic (`z_mirnov`, `z_saddle`, `z_ece`, `z_bolo`, `z_ip`,
  `z_nel`), multi-diagnostic (`z_multi`), et régression logistique
  (précurseur = 1 s avant t_tq).
- Alarme = score ≥ seuil sur 3 fenêtres consécutives ;
  seuil calibré sur le train (FA ≤ 5% des tirs sains).
- Sur tir disruptif, seules les fenêtres antérieures à t_tq comptent.

![courbes](alerte_vs_fausses_alarmes.png)

## Points de fonctionnement (seuil calibré sur train)

| Détecteur | Seuil | Fausses alarmes (test) | Détection | Alerte médiane (s) | Alerte p10 (s) |
|---|---|---|---|---|---|
| z_mirnov | 3.646005132 | 7% | 100% | 2.24 | 1.41 |
| z_saddle | 6.477699749 | 0% | 100% | 1.98 | 1.29 |
| z_ece | 5.313134634 | 0% | 0% | nan | nan |
| z_bolo | 9.637082203 | 0% | 0% | nan | nan |
| z_ip | 4.297906165 | 7% | 0% | nan | nan |
| z_nel | 1.767791888 | 7% | 27% | 1.55 | 1.20 |
| z_multi | 7.019398111 | 0% | 100% | 2.17 | 1.39 |
| logistique | 0.06844305429 | 0% | 100% | 1.72 | 1.06 |

Nota : `z_mirnov` (7%), `z_ip` (7%), `z_nel` (7%) dépassent la cible de 5% hors échantillon (15 tirs sains au test : incertitude binomiale large).

## Mode verrouillé : Mirnov vs boucles à selle (phase = 2)

Fraction des fenêtres pré-TQ au-dessus du seuil, par phase ;
« phase 2 stricte » = fenêtres entièrement postérieures au verrouillage :

| Détecteur | Phase 1 (mode tournant) | Phase 2 (toutes) | Phase 2 stricte |
|---|---|---|---|
| z_mirnov | 97% | 0% | 0% |
| z_saddle | 100% | 100% | 100% |
| z_ece | 1% | 2% | 2% |
| z_bolo | 0% | 0% | 0% |
| z_ip | 0% | 0% | 0% |
| z_nel | 24% | 18% | 17% |
| z_multi | 100% | 100% | 100% |
| logistique | 3% | 85% | 86% |

Feature dominante en phase 2 stricte : `n1_amp` (100%). Durée de la phase
verrouillée : médiane 1.77 s, max 3.14 s.

Lecture : en phase verrouillée stricte, `z_mirnov` tient 0% des fenêtres (le dB/dt s'éteint avec la rotation) tandis que `z_saddle` en tient 100% : le champ radial
statique de l'îlot verrouillé continue de croître, exactement ce que
mesure un détecteur de mode verrouillé réel. C'est la réponse
instrumentale au trou aveugle magnétique : il ne se comble pas par de
la mémoire d'alarme mais par le BON capteur. Le multi-diagnostic
(`z_multi` : 100%) et la logistique (86%) en héritent.

![mode verrouillé](mode_verrouille.png)

## Reproduire

```
node scripts/generate_v6_diags.js --shots 100 --out data/run_v6_diag02 --disrupt-ratio 0.5 --seed-base 3000
python ml/features_v6.py data/run_v6_diag02
python ml/baseline_v6.py data/run_v6_diag02
python ml/evaluate_v6.py data/run_v6_diag02
```
