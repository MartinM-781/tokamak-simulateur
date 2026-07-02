# Résultats baseline — détection précoce de disruptions

Dataset : `data/run01` — test : 30 tirs disruptifs, 30 tirs sains (split PAR TIR, 70/30 stratifié, graine fixe).

## Protocole

- Features par fenêtres glissantes de 20 ms (pas 5 ms), horodatées fin de fenêtre,
  sur les seuls canaux mesurés : RMS Mirnov, fréquence par passages à zéro,
  pente de T_e, écart d'I_p à sa médiane mobile, et leurs dérivées.
- Détecteurs z-score robustes (médiane/MAD appris sur les tirs sains du train)
  mono-canal (`z_mirnov`, `z_te`, `z_ip`) et multi-canal (`z_multi` = max des z),
  plus une régression logistique (précurseur = fenêtre dans les 300 ms avant t_tq).
- Alarme = score ≥ seuil sur 3 fenêtres consécutives. Seuil calibré
  sur le train (fausses alarmes ≤ 5% des tirs sains du train).
- Sur tir disruptif, seules les fenêtres antérieures à t_tq comptent : une
  « détection » après le quench thermique est un échec (manque).

![courbes](alerte_vs_fausses_alarmes.png)

## Points de fonctionnement (seuil calibré sur train)

| Détecteur | Seuil | Fausses alarmes (test) | Détection | Alerte médiane (ms) | Alerte p10 (ms) |
|---|---|---|---|---|---|
| z_mirnov | 12.9 | 0% | 97% | 87.7 | 57.3 |
| z_te | 2.46 | 0% | 0% | nan | nan |
| z_ip | 6.16 | 0% | 0% | nan | nan |
| z_prad | 4.25 | 3% | 100% | 67.7 | 56.2 |
| z_ne | 5.74 | 0% | 0% | nan | nan |
| z_multi | 12.9 | 0% | 100% | 87.0 | 39.4 |
| logistique | 0.713 | 10% | 100% | 110.3 | 84.6 |

## Comportement pendant le mode verrouillé (phase = 2, Mirnov muette)

Le signal Mirnov est ∝ W²·Ω : au verrouillage (Ω → 0), la sonde se tait
alors que l'îlot continue de croître. Fraction des fenêtres pré-TQ
au-dessus du seuil de fonctionnement, par phase (tirs disruptifs du test) ;
« phase 2 stricte » = fenêtres entièrement postérieures au verrouillage
(sans chevauchement du transitoire d'effondrement) :

| Détecteur | Phase 1 (précurseur rotatif) | Phase 2 (toutes fenêtres) | Phase 2 stricte |
|---|---|---|---|
| z_mirnov | 70% | 50% | 20% |
| z_te | 1% | 2% | 2% |
| z_ip | 0% | 0% | 0% |
| z_prad | 57% | 100% | 100% |
| z_ne | 0% | 0% | 0% |
| z_multi | 70% | 84% | 82% |
| logistique | 99% | 100% | 99% |

Trois enseignements (feature dominante en phase 2 : `prad_mean` (45%) puis `d_rms_mirnov` (41%) ; trou de phase 2 : 45 ms en médiane, max 55 ms, juste avant le quench) :

1. **L'alerte s'acquiert en phase 1, portée par le RMS Mirnov du mode
   tournant** (70% des fenêtres au-dessus du seuil) —
   mais les bouffées ELM des tirs sains gonflent la calibration du canal
   magnétique et lui coûtent de la marge : son seuil monte, son alerte
   recule, et il peut manquer les tirs à précurseur court.
2. **Au verrouillage, le canal magnétique s'éteint** : en phase 2 stricte
   `z_mirnov` ne tient plus que 20% des fenêtres — c'est le trou
   aveugle historique de la détection purement magnétique. Ni T_e (dérive
   lente noyée dans le clapotis des dents de scie qui a calibré la MAD),
   ni I_p (plat jusqu'au quench), ni la densité ne prennent le relais
   (`z_te`, `z_ip`, `z_ne` : ~0 % partout à 5 % de fausses alarmes).
3. **C'est P_rad qui comble le trou** : l'îlot verrouillé continue de
   grossir et de rayonner, `z_prad` tient 100% des fenêtres de
   phase 2 stricte, et le multi-canal (82%) comme la logistique
   (99%) restent armés du verrouillage jusqu'au quench.

C'est la réponse à la question mono vs multi : un détecteur mono-canal
magnétique sans mémoire perd l'alarme pendant les dizaines de
millisecondes qui précèdent immédiatement le quench (il ne survivrait que
par verrouillage d'alarme, ou par une feature à mémoire du type « le RMS
s'est effondré depuis un niveau élevé », que `d_rms_mirnov` capture au
transitoire). Le multi-canal, lui, n'a pas besoin de mémoire : la
complémentarité physique des canaux (mirnov pour le mode tournant, P_rad
pour l'îlot verrouillé) assure une couverture continue.

![mode verrouillé](mode_verrouille.png)

## Reproduire

```
node scripts/generate.js --shots 200 --out data/run01 --disrupt-ratio 0.5 --seed-base 1000
python ml/features.py data/run01
python ml/baseline.py data/run01
python ml/evaluate.py data/run01
```
