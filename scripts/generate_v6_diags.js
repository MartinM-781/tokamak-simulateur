#!/usr/bin/env node
// Générateur de datasets v6 « diagnostics » (palier B) : signaux de salle de
// contrôle multi-canaux (réseau Mirnov, boucles à selle, ECE radial,
// bolométrie, interférométrie, Ip) échantillonnés à haute cadence, un
// CSV.gz par tir + manifest.json avec la géométrie des capteurs.
//
//   node scripts/generate_v6_diags.js --shots 40 --out data/run_v6_diag01 \
//        --disrupt-ratio 0.5 --seed-base 3000 --rate 10000
//
// La rotation initiale est tirée dans [1, 2.8] kHz : avec l'harmonique 3/2 à
// 1.5×f, le contenu spectral reste sous 45 % de la cadence (anti-repliement).
'use strict';
const fs = require('fs');
const path = require('path');
const zlib = require('zlib');
const M = require(path.join(__dirname, '..', 'model', 'tokamak_model_v6.js'));
const Dg = require(path.join(__dirname, '..', 'model', 'tokamak_diags_v6.js'));

const RANGES = {
  d0a_disruptif: [0.5, 2.5],
  d0a_sain: [-1, -0.2],
  f0k: [1, 2.8],      // kHz — borné par l'anti-repliement (voir en-tête)
  cwf: [0.5, 2],
  pheat: [6, 20],     // MW
  n19: [2, 4.5],
  noise: [0.5, 3],    // %
};

function usage(msg) {
  if (msg) console.error('Erreur : ' + msg + '\n');
  console.error('Usage : node scripts/generate_v6_diags.js --out data/run_v6_diagXX ' +
                '[--shots 40] [--disrupt-ratio 0.5] [--seed-base 3000] [--rate 10000]');
  process.exit(1);
}

function parseArgs(argv) {
  const args = { shots: 40, out: null, disruptRatio: 0.5, seedBase: 3000, rate: 10000 };
  for (let i = 2; i < argv.length; i += 2) {
    const k = argv[i], v = argv[i + 1];
    if (v === undefined) usage('valeur manquante pour ' + k);
    if (k === '--shots') args.shots = Number(v);
    else if (k === '--out') args.out = v;
    else if (k === '--disrupt-ratio') args.disruptRatio = Number(v);
    else if (k === '--seed-base') args.seedBase = Number(v);
    else if (k === '--rate') args.rate = Number(v);
    else usage('option inconnue : ' + k);
  }
  if (!args.out) usage('--out est obligatoire');
  if (!Number.isInteger(args.shots) || args.shots < 1) usage('--shots doit être un entier >= 1');
  if (!(args.disruptRatio >= 0 && args.disruptRatio <= 1)) usage('--disrupt-ratio doit être dans [0, 1]');
  if (!Number.isInteger(args.seedBase)) usage('--seed-base doit être un entier');
  if (!(args.rate >= 1000)) usage('--rate doit être >= 1000 Hz');
  // Anti-repliement : l'harmonique la plus rapide (3/2 à 1.5×f0) doit rester
  // sous 45 % de la cadence d'échantillonnage.
  const fMaxHz = RANGES.f0k[1] * 1e3 * 1.5;
  if (fMaxHz > 0.45 * args.rate) {
    usage('cadence trop basse : f0 max ' + RANGES.f0k[1] + ' kHz (harmonique 3/2 à ' +
          fMaxHz + ' Hz) exige rate >= ' + Math.ceil(fMaxHz / 0.45) + ' Hz');
  }
  return args;
}

function tire(rng, [lo, hi], decimales) {
  return Number((lo + rng() * (hi - lo)).toFixed(decimales));
}

// Formats par canal : t, 8 bobines (T/s), 4 selles (T, 7 déc.), 6 ECE (keV),
// bolo (MW), nel (1e19 m^-2), ip (MA) — puis vérité terrain.
function formats() {
  const f = [5];
  for (let k = 0; k < Dg.DCFG.NCOILS; k++) f.push(4);
  for (let k = 0; k < Dg.DCFG.NSADDLE; k++) f.push(7);
  for (let k = 0; k < Dg.DCFG.ECER.length; k++) f.push(4);
  f.push(2, 4, 4);
  return f;
}

function runShot(P, seed, rate) {
  const rng = M.mulberry32(seed), g = M.makeGauss(rng), S = M.newState(P, rng);
  const p = M.defP(P), D = Dg.newDiag(P, S);
  const dts = 1 / rate, fmt = formats();
  const lignes = [];
  let tSample = 0;
  while (S.t < M.C.TMAX && !(S.ended && S.t > S.tEnd + 0.1)) {
    M.stepModel(S, P, g, M.chooseDt(S));
    if (S.t >= tSample) {
      const row = Dg.sampleDiag(D, S, p, g, dts);
      const parts = row.map((v, i) => v.toFixed(fmt[i]));
      parts.push((S.w).toFixed(5), (S.w32).toFixed(5), S.K.toFixed(4),
                 (S.Om / (2 * Math.PI) / 1e3).toFixed(4), S.locked ? 1 : 0, S.phaseId);
      lignes.push(parts.join(','));
      tSample += dts;
    }
  }
  return { lignes, S };
}

function main() {
  const args = parseArgs(process.argv);
  const debut = process.hrtime.bigint();
  const meta = M.mulberry32(args.seedBase);
  const MACH = M.resolveMachine();

  const nDisruptifs = Math.round(args.shots * args.disruptRatio);
  const labels = Array.from({ length: args.shots }, (_, i) => i < nDisruptifs);
  for (let i = labels.length - 1; i > 0; i--) {
    const j = Math.floor(meta() * (i + 1));
    const tmp = labels[i]; labels[i] = labels[j]; labels[j] = tmp;
  }

  fs.mkdirSync(args.out, { recursive: true });
  const entete = Dg.channelNames().concat(['w_m', 'w32_m', 'k_chirikov', 'f_khz', 'locked', 'phase']);
  const manifest = {
    schema: 'tokamak_v6_diag_run/1',
    generateur: 'scripts/generate_v6_diags.js',
    machine: { preset: 'jet-like', R0_m: MACH.R0, a_m: MACH.A, B0_T: MACH.B0, Ip0_MA: MACH.IP0 / 1e6 },
    rate_hz: args.rate,
    colonnes: entete,
    diagnostics: {
      mirnov_thetas_deg: Array.from({ length: Dg.DCFG.NCOILS }, (_, k) => 360 * k / Dg.DCFG.NCOILS),
      mirnov_rayon_capteur_sur_a: Dg.DCFG.BWALL,
      saddle_phis_deg: Array.from({ length: Dg.DCFG.NSADDLE }, (_, k) => 360 * k / Dg.DCFG.NSADDLE),
      ece_rayons_sur_a: Dg.DCFG.ECER,
      ece_rayon_inversion_sur_a: Dg.DCFG.RINV,
      bolo_tau_s: Dg.DCFG.TAUBOLO,
      unites: { mir: 'T/s', sad: 'T', ece: 'keV', bolo: 'MW', nel: '1e19 m^-2', ip: 'MA' },
    },
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
    const { lignes, S } = runShot(P, seed, args.rate);
    const fichier = 'tir_' + seed + '_tokamak_v6d.csv.gz';
    const csv = '# tokamak_v6 diag shot (SI, JET-like)\n' +
      '# seed=' + seed + ' d0a=' + P.d0a + ' f0_khz=' + P.f0k + ' cwf=' + P.cwf +
      ' pheat_mw=' + P.pheat + ' n19=' + P.n19 + ' noise_pct=' + P.noise + '\n' +
      '# t_lock_s=' + (S.tLock >= 0 ? S.tLock.toFixed(4) : 'NA') +
      ' t_tq_s=' + (S.tTQ >= 0 ? S.tTQ.toFixed(4) : 'NA') +
      ' t_cq_s=' + (S.tCQ >= 0 ? S.tCQ.toFixed(4) : 'NA') +
      ' t_onset_s=' + (S.tOnset >= 0 ? S.tOnset.toFixed(4) : 'NA') + '\n' +
      entete.join(',') + '\n' + lignes.join('\n') + '\n';
    fs.writeFileSync(path.join(args.out, fichier), zlib.gzipSync(Buffer.from(csv), { level: 6 }));

    const disrupte = S.tTQ >= 0;
    if (disrupte !== disruptif) incoherences++;
    manifest.tirs.push({
      fichier, seed,
      label: disruptif ? 'disruptif' : 'sain',
      disrupte,
      params: P,
      t_onset_s: S.tOnset >= 0 ? S.tOnset : null,
      t_lock_s: S.tLock >= 0 ? S.tLock : null,
      t_tq_s: S.tTQ >= 0 ? S.tTQ : null,
      t_cq_s: S.tCQ >= 0 ? S.tCQ : null,
      t_end_s: S.tEnd >= 0 ? S.tEnd : null,
      w_at_lock_m: S.wAtLock >= 0 ? S.wAtLock : null,
      n_lignes_csv: lignes.length,
    });
    if ((i + 1) % 10 === 0) console.log((i + 1) + '/' + args.shots + ' tirs générés…');
  }

  fs.writeFileSync(path.join(args.out, 'manifest.json'), JSON.stringify(manifest, null, 2) + '\n');
  const dureeS = Number(process.hrtime.bigint() - debut) / 1e9;
  console.log('Terminé : ' + args.shots + ' tirs (' + nDisruptifs + ' disruptifs visés), ' +
              args.rate + ' Hz, dans ' + args.out + ' en ' + dureeS.toFixed(1) + ' s.');
  if (incoherences > 0) {
    console.warn('Attention : ' + incoherences + ' tir(s) dont l\'issue réelle contredit le label visé.');
  }
}

main();
