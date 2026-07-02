# Tokamak Simulateur

Suite pédagogique de 5 pages HTML autonomes simulant la physique des tokamaks
(v1 lignes de champ et facteur q, v2 orbites de particules, v3 équilibre MHD,
v4 îlots magnétiques, v5 dynamique de disruption complète), plus un générateur
de datasets synthétiques labellisés destinés à l'entraînement ML de prédiction
de disruptions par détection d'anomalies sur séries temporelles.

## Arborescence

- `web/` — les pages HTML autonomes. S'ouvrent directement en `file://`, sans serveur.
- `model/tokamak_model.js` — modèle physique pur de la v5 (sans DOM), double export.
- `scripts/` — outillage Node (génération de datasets en CLI).
- `tests/` — tests de régression physique (`node:test`).
- `ml/` — baseline Python (features, détecteurs, évaluation).
- `data/` — datasets générés (gitignoré, régénérable à l'identique par seed).

## Conventions — côté web

- Vanilla JS uniquement. Aucun framework, aucun bundler, aucune étape de build,
  aucune dépendance npm côté web.
- Les pages doivent TOUJOURS rester ouvrables directement en `file://` :
  pas de modules ES, pas de `fetch` de ressources locales, chemins relatifs
  uniquement pour les scripts partagés.
- Design system commun intouchable : palette définie dans `:root` (`--bg`,
  `--panel`, `--screen`, `--line`, `--text`, `--sec`, `--mut`, `--plasma`,
  `--reson`, `--warn`, `--steel`, `--cuivre`), typographies IBM Plex Mono
  (texte/UI) et Space Grotesk (titres). Seule ressource externe autorisée :
  Google Fonts.
- La physique des v1 à v4 est validée : ne pas y toucher.

## Conventions — modèle et physique

- Unités arbitraires assumées (temps en ms) : c'est un générateur pédagogique
  de précurseurs plausibles, pas un code de production.
- Le modèle v5 vit dans `model/tokamak_model.js`, entre les marqueurs
  `/*MODEL-BEGIN*/` et `/*MODEL-END*/` : constantes `MP`, `mulberry32`,
  `makeGauss`, `newState`, `stepModel`. Double export : `window.TokamakModel`
  en navigateur ET `module.exports` en Node — pas de module ES.
- `dt = 0.05 ms` partout (page web, tests, générateur).
- Le bruit est purement observationnel (AR(1) sur les mesures) : la dynamique
  d'état est déterministe à paramètres donnés ; toute la stochasticité est
  reproductible par seed (`mulberry32`).
- Toute modification de la physique doit soit passer `npm test`, soit mettre
  à jour les tests avec justification explicite dans le message de commit.
  `npm test` est obligatoire avant tout commit touchant la physique.

## Outillage

- Node ≥ 18 (aucune dépendance npm), Python ≥ 3.10 (venv dans `ml/.venv`).
- `npm test` — tests de régression physique.
- `node scripts/generate.js --shots N --out data/runXX --disrupt-ratio R --seed-base S`
  — génération de datasets (un CSV par tir + `manifest.json`).

## Commits

- Atomiques, messages en français à l'impératif.
