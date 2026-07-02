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
| z_mirnov | 5.91 | 3% | 100% | 114.3 | 85.8 |
| z_te | 2.47 | 0% | 0% | nan | nan |
| z_ip | 5.08 | 10% | 0% | nan | nan |
| z_multi | 5.9 | 7% | 100% | 114.3 | 85.8 |
| logistique | 0.617 | 3% | 100% | 114.4 | 85.8 |

## Comportement pendant le mode verrouillé (phase = 2, Mirnov muette)

Le signal Mirnov est ∝ W²·Ω : au verrouillage (Ω → 0), la sonde se tait
alors que l'îlot continue de croître. Fraction des fenêtres pré-TQ
au-dessus du seuil de fonctionnement, par phase (tirs disruptifs du test) ;
« phase 2 stricte » = fenêtres entièrement postérieures au verrouillage
(sans chevauchement du transitoire d'effondrement) :

| Détecteur | Phase 1 (précurseur rotatif) | Phase 2 (toutes fenêtres) | Phase 2 stricte |
|---|---|---|---|
| z_mirnov | 100% | 57% | 23% |
| z_te | 2% | 0% | 0% |
| z_ip | 0% | 0% | 0% |
| z_multi | 100% | 57% | 23% |
| logistique | 100% | 60% | 29% |

Trois enseignements (z_multi médian : 25 en phase 1, 19 en phase 2, 2 en phase 2 stricte) :

1. **Le 100 % de détection ne doit rien à la phase verrouillée** : l'alarme
   se déclenche en phase 1, portée par le RMS Mirnov du mode tournant, et
   le temps d'alerte est déjà acquis quand le mode se verrouille.
2. **En phase 2 stricte, tous les détecteurs fenêtrés deviennent quasi
   aveugles** (23% de dépassements pour le multi-canal,
   z médian ≈ 2) : le Mirnov retombe au niveau du bruit, la
   lente dérive de T_e est noyée dans le clapotis des dents de scie qui a
   calibré la MAD, et I_p ne bouge qu'au quench (`z_te` et `z_ip` ne
   détectent d'ailleurs RIEN à 5 % de fausses alarmes, phase 1 comprise).
   Le multi-canal n'est pas sauvé par te/ip — il retombe avec le mirnov.
   Ce trou aveugle dure 45 ms en médiane (max 55 ms), juste avant le quench.
3. **Le seul signal propre à la phase 2 est le transitoire d'extinction
   lui-même** : sur ces fenêtres, la feature dominante est `d_rms_mirnov` (51%) puis `rms_mirnov` (20%) —
   la disparition brutale du signal précurseur est elle-même un précurseur.

Conséquence pour la détection d'anomalies : un détecteur par fenêtre sans
mémoire perdrait l'alarme pendant les dizaines de millisecondes qui
précèdent immédiatement le quench. Il faut soit une alarme à verrouillage
(latching, utilisé ici), soit des features à mémoire longue capables
d'encoder « le RMS s'est effondré depuis un niveau élevé » — c'est
exactement ce que la dérivée `d_rms_mirnov` capture au moment de la
transition, et ce qu'un modèle séquentiel (RNN/TCN) apprendrait de lui-même.

![mode verrouillé](mode_verrouille.png)

## Reproduire

```
node scripts/generate.js --shots 200 --out data/run01 --disrupt-ratio 0.5 --seed-base 1000
python ml/features.py data/run01
python ml/baseline.py data/run01
python ml/evaluate.py data/run01
```
