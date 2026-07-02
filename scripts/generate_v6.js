#!/usr/bin/env node
// Générateur de datasets v6 (unités SI, machine JET-like) — même contrat que
// generate.js : reproductible de bout en bout, un CSV par tir + manifest.json
// sans horodatage.
//
//   node scripts/generate_v6.js --shots 100 --out data/run_v6_01 --disrupt-ratio 0.5 --seed-base 2000
//
// Échantillonnage de sortie : 1 kHz (le pas d'intégration est adaptatif,
// 100 µs hors quench et 10 µs pendant — voir chooseDt du modèle).
'use strict';
const fs = require('fs');
const path = require('path');
const M = require(path.join(__dirname, '..', 'model', 'tokamak_model_v6.js'));

const SAMPLE_S = 1e-3;

const RANGES = {
  d0a_disruptif: [0.5, 2.5],
  d0a_sain: [-1, -0.2],
  f0k: [3, 7],        // kHz
  cwf: [0.5, 2],      // facteur de freinage paroi
  pheat: [6, 20],     // MW
  n19: [2, 4.5],      // 1e19 m^-3
  noise: [0.5, 3],    // %
};

function usage(msg) {
  if (msg) console.error('Erreur : ' + msg + '\n');
  console.error('Usage : node scripts/generate_v6.js --out data/run_v6_XX [--shots 100] [--disrupt-ratio 0.5] [--seed-base 2000]');
  process.exit(1);
}

function parseArgs(argv) {
  const args = { shots: 100, out: null, disruptRatio: 0.5, seedBase: 2000 };
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

function tire(rng, [lo, hi], decimales) {
  return Number((lo + rng() * (hi - lo)).toFixed(decimales));
}

// Boucle de tir à pas adaptatif, échantillonnée à 1 kHz.
function runShot(P, seed) {
  const rng = M.mulberry32(seed), g = M.makeGauss(rng), S = M.newState(P, rng);
  const A = { t: [], m: [], te: [], ip: [], pr: [], ne: [], w: [], w32: [], k: [], f: [], lk: [], ph: [] };
  let tSample = 0;
  while (S.t < M.C.TMAX && !(S.ended && S.t > S.tEnd + 0.1)) {
    const meas = M.stepModel(S, P, g, M.chooseDt(S));
    if (S.t >= tSample) {
      A.t.push(S.t); A.m.push(meas.mir); A.te.push(meas.te); A.ip.push(meas.ip);
      A.pr.push(meas.prad); A.ne.push(meas.ne);
      A.w.push(S.w); A.w32.push(S.w32); A.k.push(S.K);
      A.f.push(S.Om / (2 * Math.PI) / 1e3); A.lk.push(S.locked ? 1 : 0); A.ph.push(S.phaseId);
      tSample += SAMPLE_S;
    }
  }
  return { A, S };
}

function buildCsv(seed, P, S, A) {
  const L = [];
  L.push('# tokamak_v6 synthetic shot (SI, JET-like)');
  L.push('# seed=' + seed + ' d0a=' + P.d0a + ' f0_khz=' + P.f0k + ' cwf=' + P.cwf +
         ' pheat_mw=' + P.pheat + ' n19=' + P.n19 + ' noise_pct=' + P.noise);
  L.push('# t_lock_s=' + (S.tLock >= 0 ? S.tLock.toFixed(4) : 'NA') +
         ' t_tq_s=' + (S.tTQ >= 0 ? S.tTQ.toFixed(4) : 'NA') +
         ' t_cq_s=' + (S.tCQ >= 0 ? S.tCQ.toFixed(4) : 'NA') +
         ' t_onset_s=' + (S.tOnset >= 0 ? S.tOnset.toFixed(4) : 'NA'));
  L.push('t_s,mirnov_Tps,te0_keV,ip_MA,prad_MW,ne_1e19,w_m,w32_m,k_chirikov,f_khz,locked,phase');
  for (let i = 0; i < A.t.length; i++) {
    L.push(A.t[i].toFixed(4) + ',' + A.m[i].toFixed(4) + ',' + A.te[i].toFixed(4) + ',' +
           A.ip[i].toFixed(4) + ',' + A.pr[i].toFixed(2) + ',' + A.ne[i].toFixed(4) + ',' +
           A.w[i].toFixed(5) + ',' + A.w32[i].toFixed(5) + ',' + A.k[i].toFixed(4) + ',' +
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
  for (let i = labels.length - 1; i > 0; i--) {
    const j = Math.floor(meta() * (i + 1));
    const tmp = labels[i]; labels[i] = labels[j]; labels[j] = tmp;
  }

  fs.mkdirSync(args.out, { recursive: true });
  const MACH = M.resolveMachine();
  const manifest = {
    schema: 'tokamak_v6_run/1',
    generateur: 'scripts/generate_v6.js',
    machine: { preset: 'jet-like', R0_m: MACH.R0, a_m: MACH.A, B0_T: MACH.B0, Ip0_MA: MACH.IP0 / 1e6 },
    sample_s: SAMPLE_S,
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
      d0a: tire(meta, disruptif ? RANGES.d0a_disruptif : RANGES.d0a_sain, 3),
      f0k: tire(meta, RANGES.f0k, 3),
      cwf: tire(meta, RANGES.cwf, 3),
      pheat: tire(meta, RANGES.pheat, 2),
      n19: tire(meta, RANGES.n19, 3),
      noise: tire(meta, RANGES.noise, 3),
    };
    const seed = args.seedBase + i;
    const { A, S } = runShot(P, seed);
    const fichier = 'tir_' + seed + '_tokamak_v6.csv';
    fs.writeFileSync(path.join(args.out, fichier), buildCsv(seed, P, S, A) + '\n');

    const disrupte = S.tTQ >= 0;
    if (disrupte !== disruptif) incoherences++;
    manifest.tirs.push({
      fichier,
      seed,
      label: disruptif ? 'disruptif' : 'sain',
      disrupte,
      params: P,
      t_onset_s: S.tOnset >= 0 ? S.tOnset : null,
      t_lock_s: S.tLock >= 0 ? S.tLock : null,
      t_tq_s: S.tTQ >= 0 ? S.tTQ : null,
      t_cq_s: S.tCQ >= 0 ? S.tCQ : null,
      t_end_s: S.tEnd >= 0 ? S.tEnd : null,
      w_at_lock_m: S.wAtLock >= 0 ? S.wAtLock : null,
      n_lignes_csv: A.t.length,
    });
    if ((i + 1) % 25 === 0) console.log((i + 1) + '/' + args.shots + ' tirs générés…');
  }

  fs.writeFileSync(path.join(args.out, 'manifest.json'), JSON.stringify(manifest, null, 2) + '\n');
  const dureeS = Number(process.hrtime.bigint() - debut) / 1e9;
  console.log('Terminé : ' + args.shots + ' tirs (' + nDisruptifs + ' disruptifs visés) dans ' +
              args.out + ' en ' + dureeS.toFixed(1) + ' s.');
  if (incoherences > 0) {
    console.warn('Attention : ' + incoherences + ' tir(s) dont l\'issue réelle contredit le label visé.');
  }
}

main();
