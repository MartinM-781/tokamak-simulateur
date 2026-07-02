// Tests de la couche diagnostics v6 (palier B) — les signatures qu'un
// physicien attend d'une vraie salle de contrôle :
//   structure poloïdale m=2 dans le réseau Mirnov, boucles à selle qui
//   continuent de voir le mode verrouillé quand les Mirnov se taisent,
//   dents de scie inversées autour du rayon d'inversion sur l'ECE,
//   aplatissement du profil au passage de l'îlot, bolomètre retardé.
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const { execFileSync } = require('node:child_process');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const zlib = require('node:zlib');
const M = require('../model/tokamak_model_v6.js');
const Dg = require('../model/tokamak_diags_v6.js');

const GEN = path.join(__dirname, '..', 'scripts', 'generate_v6_diags.js');
// Indices de colonnes : t=0, bobines 1-8, selles 9-12, ECE 13-18, bolo 19,
// nel 20, ip 21 (voir channelNames).
const COIL0 = 1, SAD0 = 9, ECE0 = 13, BOLO = 19;
const RATE = 10000;

function runDiag(P, seed, tStop) {
  const rng = M.mulberry32(seed), g = M.makeGauss(rng), S = M.newState(P, rng);
  const p = M.defP(P), D = Dg.newDiag(P, S);
  const dts = 1 / RATE;
  const rows = [], te0 = [], pradMW = [];
  let tSample = 0;
  while (S.t < tStop && !(S.ended && S.t > S.tEnd + 0.05)) {
    M.stepModel(S, P, g, M.chooseDt(S));
    if (S.t >= tSample) {
      rows.push(Dg.sampleDiag(D, S, p, g, dts));
      te0.push(S.Te0);
      pradMW.push((S.Prad + S.prFlash) / 1e6);
      tSample += dts;
    }
  }
  return { rows, te0, pradMW, S };
}

const cache = new Map();
function shot(P, seed, tStop) {
  const k = JSON.stringify([P, seed, tStop]);
  if (!cache.has(k)) cache.set(k, runDiag(P, seed, tStop));
  return cache.get(k);
}
const CLASSIQUE = { d0a: 1, f0k: 2 };

// Amplitude du nombre de mode poloïdal m par DFT sur le réseau de bobines.
function ampM(row, m) {
  let re = 0, im = 0;
  for (let j = 0; j < Dg.DCFG.NCOILS; j++) {
    const th = 2 * Math.PI * j / Dg.DCFG.NCOILS;
    re += row[COIL0 + j] * Math.cos(m * th);
    im += row[COIL0 + j] * Math.sin(m * th);
  }
  return Math.hypot(re, im) / Dg.DCFG.NCOILS;
}
// Amplitude n=1 reconstruite des 4 boucles à selle (différences opposées).
function ampN1(row) {
  return Math.hypot(row[SAD0] - row[SAD0 + 2], row[SAD0 + 1] - row[SAD0 + 3]) / 2;
}
const rms = a => Math.sqrt(a.reduce((s, x) => s + x * x, 0) / a.length);
const moy = a => a.reduce((s, x) => s + x, 0) / a.length;

test('réseau Mirnov : la structure poloïdale m=2 domine pendant le précurseur', () => {
  const o = shot(CLASSIQUE, 42, 3.6);
  const { S, rows } = o;
  assert.ok(S.tLock > 0, 'le tir de référence doit se verrouiller');
  const i0 = Math.floor((S.tLock - 0.4) * RATE), i1 = Math.floor((S.tLock - 0.1) * RATE);
  let m1 = 0, m2 = 0, m3 = 0, n = 0;
  for (let i = i0; i < i1; i += 7) { m1 += ampM(rows[i], 1); m2 += ampM(rows[i], 2); m3 += ampM(rows[i], 3); n++; }
  m1 /= n; m2 /= n; m3 /= n;
  assert.ok(m2 > 2 * m3, `m2=${m2} doit dominer m3=${m3}`);
  assert.ok(m2 > 4 * m1, `m2=${m2} doit dominer m1=${m1}`);
});

test('boucles à selle : le mode verrouillé reste visible quand les Mirnov se taisent', () => {
  const o = shot(CLASSIQUE, 42, 3.6);
  const { S, rows } = o;
  const iTot = Math.floor(0.2 * RATE), iTot1 = Math.floor(0.5 * RATE);
  const iLk = Math.floor((S.tLock + 0.05) * RATE), iLk1 = Math.floor((S.tTQ - 0.01) * RATE);
  const iPh1 = Math.floor((S.tLock - 0.4) * RATE), iPh11 = Math.floor((S.tLock - 0.05) * RATE);
  const n1Tot = [], n1Lock = [], coilPh1 = [], coilLock = [];
  for (let i = iTot; i < iTot1; i++) n1Tot.push(ampN1(rows[i]));
  for (let i = iLk; i < iLk1; i++) { n1Lock.push(ampN1(rows[i])); coilLock.push(rows[i][COIL0]); }
  for (let i = iPh1; i < iPh11; i++) coilPh1.push(rows[i][COIL0]);
  // La selle voit l'îlot verrouillé grossir…
  assert.ok(moy(n1Lock) > 5 * moy(n1Tot),
    `n1 verrouillé=${moy(n1Lock)} vs début de tir=${moy(n1Tot)}`);
  // …pendant que la bobine dB/dt retombe vers son plancher de bruit.
  assert.ok(rms(coilLock) < 0.3 * rms(coilPh1),
    `rms bobine verrouillé=${rms(coilLock)} vs précurseur=${rms(coilPh1)}`);
});

test('ECE : dents de scie inversées — le cœur chute, l’extérieur reçoit le pulse', () => {
  const o = shot({ d0a: -0.5, f0k: 2 }, 42, 1.8);
  const { rows, te0 } = o;
  const crashs = [];
  for (let i = 1; i < te0.length; i++) if (te0[i] < te0[i - 1] * 0.93) crashs.push(i);
  assert.ok(crashs.length >= 5, `crashs détectés=${crashs.length}`);
  let ok = 0, n = 0;
  for (const i of crashs) {
    if (i < 30 || i + 30 >= rows.length) continue;
    n++;
    const avant = k => moy(rows.slice(i - 25, i - 2).map(r => r[ECE0 + k]));
    const apres = k => moy(rows.slice(i + 2, i + 25).map(r => r[ECE0 + k]));
    const coeur = apres(1) - avant(1);      // ece_r025 : doit chuter
    const bord = apres(5) - avant(5);       // ece_r085 : doit monter (pulse)
    if (coeur < 0 && bord > 0) ok++;
  }
  assert.ok(ok >= 0.7 * n, `inversion vue sur ${ok}/${n} crashs`);
});

test('ECE : l’îlot aplatit le profil au voisinage de q=2 (test unitaire, sans bruit)', () => {
  const P = { noise: 0, elm: 0 };
  const p = M.defP(P);
  const rng = M.mulberry32(7), g = M.makeGauss(rng);
  const S = M.newState(P, rng);
  S.Te0 = 5; S.w = 0.12;                       // îlot de 12 cm posé sur q=2
  const D = Dg.newDiag(P, S);
  let row;
  for (let i = 0; i < 6000; i++) row = Dg.sampleDiag(D, S, p, g, 1e-4);
  const attenduPlat = 5 * Dg.profil(0.70);     // valeur sans îlot
  const attenduLoin = 5 * Dg.profil(0.25);
  assert.ok(row[ECE0 + 4] < attenduPlat * 0.995 && row[ECE0 + 4] > attenduPlat * 0.9,
    `ece_r070=${row[ECE0 + 4]} vs profil nu=${attenduPlat}`);
  assert.ok(Math.abs(row[ECE0 + 1] - attenduLoin) < attenduLoin * 0.002,
    `ece_r025=${row[ECE0 + 1]} doit rester au profil nu=${attenduLoin}`);
});

test('bolomètre : le flash du quench est lissé et retardé par le temps de réponse', () => {
  const o = shot(CLASSIQUE, 42, 3.6);
  const { rows, pradMW, S } = o;
  const flashMax = Math.max(...pradMW);
  let boloMax = 0, tBoloMax = 0;
  for (const r of rows) if (r[BOLO] > boloMax) { boloMax = r[BOLO]; tBoloMax = r[0]; }
  assert.ok(flashMax > 1000, `flash instantané=${flashMax} MW (attendu échelle GW)`);
  assert.ok(boloMax < flashMax / 2, `bolo max=${boloMax} doit être lissé (< flash/2)`);
  assert.ok(boloMax >= 50 && boloMax <= 2000, `bolo max=${boloMax} MW`);
  assert.ok(tBoloMax >= S.tTQ, 'le pic bolométrique doit suivre le TQ');
});

test('diagnostics : déterminisme complet par seed', () => {
  const a = runDiag({ d0a: -0.5, f0k: 2 }, 11, 0.4);
  const b = runDiag({ d0a: -0.5, f0k: 2 }, 11, 0.4);
  assert.deepEqual(a.rows, b.rows);
});

test('générateur diag : manifest avec géométrie, CSV.gz relisible et déterministe', () => {
  const d1 = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_diag_a_'));
  const d2 = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_diag_b_'));
  try {
    const args = ['--shots', '2', '--disrupt-ratio', '0.5', '--seed-base', '555'];
    execFileSync(process.execPath, [GEN, '--out', d1, ...args], { encoding: 'utf8' });
    execFileSync(process.execPath, [GEN, '--out', d2, ...args], { encoding: 'utf8' });
    const man = JSON.parse(fs.readFileSync(path.join(d1, 'manifest.json'), 'utf8'));
    assert.equal(man.schema, 'tokamak_v6_diag_run/1');
    assert.equal(man.rate_hz, 10000);
    assert.equal(man.diagnostics.mirnov_thetas_deg.length, 8);
    assert.equal(man.diagnostics.ece_rayons_sur_a.length, 6);
    assert.equal(man.colonnes.length, 22 + 6);
    for (const tir of man.tirs) {
      const brut = zlib.gunzipSync(fs.readFileSync(path.join(d1, tir.fichier))).toString('utf8');
      const lignes = brut.trim().split('\n');
      assert.equal(lignes[3], man.colonnes.join(','));
      assert.equal(lignes.length - 4, tir.n_lignes_csv);
      assert.equal(lignes[4].split(',').length, man.colonnes.length);
    }
    for (const f of fs.readdirSync(d1).sort()) {
      assert.ok(fs.readFileSync(path.join(d1, f)).equals(fs.readFileSync(path.join(d2, f))),
        `divergence sur ${f}`);
    }
  } finally {
    fs.rmSync(d1, { recursive: true, force: true });
    fs.rmSync(d2, { recursive: true, force: true });
  }
});

test('générateur diag : refuse une cadence sous la limite anti-repliement', () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_diag_c_'));
  try {
    assert.throws(
      () => execFileSync(process.execPath, [GEN, '--out', dir, '--shots', '1', '--rate', '2000'],
        { encoding: 'utf8', stdio: 'pipe' }),
      /cadence trop basse/
    );
  } finally {
    fs.rmSync(dir, { recursive: true, force: true });
  }
});
