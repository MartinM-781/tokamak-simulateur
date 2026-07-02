#!/usr/bin/env node
// Générateur de datasets synthétiques de disruptions (CLI, Node >= 18, zéro dépendance).
//
//   node scripts/generate.js --shots 200 --out data/run01 --disrupt-ratio 0.5 --seed-base 1000
//
// Reproductible de bout en bout : un méta-RNG mulberry32(seed-base) tire les
// paramètres et mélange les labels ; le tir i utilise la graine seed-base + i.
// Le manifest ne contient volontairement aucun horodatage : mêmes arguments
// ⇒ sorties identiques octet pour octet.
'use strict';
const fs = require('fs');
const path = require('path');
const M = require(path.join(__dirname, '..', 'model', 'tokamak_model.js'));

const DT = 0.05;      // ms — partout (web, tests, générateur)
const CSV_STRIDE = 2; // 1 ligne sur 2 pas, comme l'export web (échantillonnage 0.1 ms)

const RANGES = {
  d0_disruptif: [0.5, 2.5],
  d0_sain: [-1, -0.2],
  f0: [2, 5],
  cb: [10, 35],
  noise: [0.5, 3],
};

function usage(msg) {
  if (msg) console.error('Erreur : ' + msg + '\n');
  console.error('Usage : node scripts/generate.js --out data/runXX [--shots 200] [--disrupt-ratio 0.5] [--seed-base 1000]');
  process.exit(1);
}

function parseArgs(argv) {
  const args = { shots: 200, out: null, disruptRatio: 0.5, seedBase: 1000 };
  for (let i = 2; i < argv.length; i += 2) {
    const k = argv[i], v = argv[i + 1];
    if (v === undefined) usage('valeur manquante pour ' + k);
    if (k === '--shots') args.shots = Number(v);
    else if (k === '--out') args.out = v;
    else if (k === '--disrupt-ratio') args.disruptRatio = Number(v);
    else if (k === '--seed-base') args.seedBase = Number(v);
    else usage('option inconnue : ' + k);
  }
  if (!args.out) usage('--out est obligatoire');
  if (!Number.isInteger(args.shots) || args.shots < 1) usage('--shots doit être un entier >= 1');
  if (!(args.disruptRatio >= 0 && args.disruptRatio <= 1)) usage('--disrupt-ratio doit être dans [0, 1]');
  if (!Number.isInteger(args.seedBase)) usage('--seed-base doit être un entier');
  return args;
}

// Tirage uniforme arrondi (l'arrondi fait partie du protocole : la valeur
// arrondie est celle utilisée par la simulation ET écrite dans le CSV).
function tire(rng, [lo, hi], decimales) {
  return Number((lo + rng() * (hi - lo)).toFixed(decimales));
}

// Boucle de tir identique à la page web : arrêt à TMAX ou 60 ms après tEnd,
// échantillon poussé après stepModel (S.t déjà incrémenté).
function runShot(P, seed) {
  const rng = M.mulberry32(seed), g = M.makeGauss(rng), S = M.newState(P);
  const A = { t: [], m: [], te: [], ip: [], w1: [], w2: [], k: [], f: [], lk: [], ph: [] };
  while (S.t < M.MP.TMAX && !(S.ended && S.t > S.tEnd + 60)) {
    const meas = M.stepModel(S, P, g, DT);
    A.t.push(S.t); A.m.push(meas.mir); A.te.push(meas.te); A.ip.push(meas.ip);
    A.w1.push(S.W21); A.w2.push(S.W32); A.k.push(S.K);
    A.f.push(S.Om / (2 * Math.PI)); A.lk.push(S.locked ? 1 : 0); A.ph.push(S.phaseId);
  }
  return { A, S };
}

// Format strictement identique au buildCsv() de web/tokamak_v5.html.
function buildCsv(seed, P, S, A) {
  const L = [];
  L.push('# tokamak_v5 synthetic shot');
  L.push('# seed=' + seed + ' delta0=' + P.d0 + ' f0_khz=' + P.f0 + ' cb=' + P.cb + ' noise_pct=' + P.noise);
  L.push('# t_lock_ms=' + (S.tLock >= 0 ? S.tLock.toFixed(2) : 'NA') +
         ' t_tq_ms=' + (S.tTQ >= 0 ? S.tTQ.toFixed(2) : 'NA') +
         ' t_cq_ms=' + (S.tCQ >= 0 ? S.tCQ.toFixed(2) : 'NA'));
  L.push('t_ms,mirnov,te,ip,w21,w32,k_chirikov,f_khz,locked,phase');
  for (let i = 0; i < A.t.length; i += CSV_STRIDE) {
    L.push(A.t[i].toFixed(2) + ',' + A.m[i].toFixed(5) + ',' + A.te[i].toFixed(5) + ',' + A.ip[i].toFixed(5) + ',' +
           A.w1[i].toFixed(5) + ',' + A.w2[i].toFixed(5) + ',' + A.k[i].toFixed(4) + ',' +
           A.f[i].toFixed(4) + ',' + A.lk[i] + ',' + A.ph[i]);
  }
  return L.join('\n');
}

function main() {
  const args = parseArgs(process.argv);
  const debut = process.hrtime.bigint();
  const meta = M.mulberry32(args.seedBase);

  const nDisruptifs = Math.round(args.shots * args.disruptRatio);
  const labels = Array.from({ length: args.shots }, (_, i) => i < nDisruptifs);
  for (let i = labels.length - 1; i > 0; i--) { // Fisher-Yates via méta-RNG
    const j = Math.floor(meta() * (i + 1));
    const tmp = labels[i]; labels[i] = labels[j]; labels[j] = tmp;
  }

  fs.mkdirSync(args.out, { recursive: true });
  const manifest = {
    schema: 'tokamak_v5_run/1',
    generateur: 'scripts/generate.js',
    dt_ms: DT,
    csv_stride: CSV_STRIDE,
    params: {
      shots: args.shots,
      disrupt_ratio: args.disruptRatio,
      seed_base: args.seedBase,
      ranges: RANGES,
    },
    tirs: [],
  };

  let incoherences = 0;
  for (let i = 0; i < args.shots; i++) {
    const disruptif = labels[i];
    const P = {
      d0: tire(meta, disruptif ? RANGES.d0_disruptif : RANGES.d0_sain, 3),
      f0: tire(meta, RANGES.f0, 3),
      cb: tire(meta, RANGES.cb, 2),
      noise: tire(meta, RANGES.noise, 3),
    };
    const seed = args.seedBase + i;
    const { A, S } = runShot(P, seed);
    const fichier = 'tir_' + seed + '_tokamak_v5.csv';
    fs.writeFileSync(path.join(args.out, fichier), buildCsv(seed, P, S, A) + '\n');

    const disrupte = S.tTQ >= 0;
    if (disrupte !== disruptif) incoherences++;
    manifest.tirs.push({
      fichier,
      seed,
      label: disruptif ? 'disruptif' : 'sain',
      disrupte,
      params: P,
      t_lock_ms: S.tLock >= 0 ? S.tLock : null,
      t_tq_ms: S.tTQ >= 0 ? S.tTQ : null,
      t_cq_ms: S.tCQ >= 0 ? S.tCQ : null,
      t_end_ms: S.tEnd >= 0 ? S.tEnd : null,
      n_pas: A.t.length,
      n_lignes_csv: Math.ceil(A.t.length / CSV_STRIDE),
    });
    if ((i + 1) % 50 === 0) console.log((i + 1) + '/' + args.shots + ' tirs générés…');
  }

  fs.writeFileSync(path.join(args.out, 'manifest.json'), JSON.stringify(manifest, null, 2) + '\n');
  const dureeS = Number(process.hrtime.bigint() - debut) / 1e9;
  console.log('Terminé : ' + args.shots + ' tirs (' + nDisruptifs + ' disruptifs visés) dans ' +
              args.out + ' en ' + dureeS.toFixed(1) + ' s.');
  if (incoherences > 0) {
    console.warn('Attention : ' + incoherences + ' tir(s) dont l\'issue réelle contredit le label visé (voir manifest, champs label/disrupte).');
  }
}

main();
