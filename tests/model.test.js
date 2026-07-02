// Tests de régression physique du modèle v5 (model/tokamak_model.js).
// Invariants validés par exécution de référence le 2026-07-02 (dt = 0.05 ms) :
//   classique seed 42 : tLock=144.40, tTQ=187.60, tCQ=188.85, tEnd=213.45 ms,
//   pic ip mesuré=1.0656 ; sain seed 42 : min S.Te=0.8314, 104 dents de scie,
//   std bruit te=0.01514 (théorie AR(1) : 0.0150).
// Toute modification physique doit soit passer ces tests, soit les mettre à
// jour avec justification explicite dans le message de commit.
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const M = require('../model/tokamak_model.js');

const DT = 0.05;
const CLASSIQUE = { d0: 1, f0: 3, cb: 20, noise: 1.5 };
const SAIN = { d0: -0.5, f0: 3, cb: 20, noise: 1.5 };

// Reproduit la boucle de tir de la page web : arrêt à TMAX ou 60 ms après tEnd.
function runShot(P, seed) {
  const rng = M.mulberry32(seed);
  const g = M.makeGauss(rng);
  const S = M.newState(P);
  const out = { t: [], mir: [], te: [], ip: [], prad: [], ne: [], teEtat: [], dents: 0 };
  let tePrec = S.Te;
  while (S.t < M.MP.TMAX && !(S.ended && S.t > S.tEnd + 60)) {
    const meas = M.stepModel(S, P, g, DT);
    out.t.push(S.t);
    out.mir.push(meas.mir);
    out.te.push(meas.te);
    out.ip.push(meas.ip);
    out.prad.push(meas.prad);
    out.ne.push(meas.ne);
    out.teEtat.push(S.Te);
    // Dent de scie : chute de S.Te de plus de 5 % en un seul pas (le crash
    // est instantané : Te *= 0.92 ; hors crash, |dTe| par pas << 1 %).
    if (S.Te < tePrec * 0.95) out.dents++;
    tePrec = S.Te;
  }
  out.S = S;
  return out;
}

test('classique (seed 42) : ordre strict des événements tLock < tTQ < tCQ < tEnd', () => {
  const { S } = runShot(CLASSIQUE, 42);
  assert.ok(S.tLock >= 0, `tLock=${S.tLock} : le mode doit se verrouiller`);
  assert.ok(S.tLock < S.tTQ, `tLock=${S.tLock} doit précéder tTQ=${S.tTQ}`);
  assert.ok(S.tTQ < S.tCQ, `tTQ=${S.tTQ} doit précéder tCQ=${S.tCQ}`);
  assert.ok(S.tCQ < S.tEnd, `tCQ=${S.tCQ} doit précéder tEnd=${S.tEnd}`);
});

test('classique (seed 42) : tLock dans [100, 200] ms', () => {
  const { S } = runShot(CLASSIQUE, 42);
  assert.ok(S.tLock >= 100 && S.tLock <= 200, `tLock=${S.tLock}`);
});

test('classique (seed 42) : pic d’Ip mesuré dans [1.03, 1.12]', () => {
  const { ip } = runShot(CLASSIQUE, 42);
  const pic = Math.max(...ip);
  assert.ok(pic >= 1.03 && pic <= 1.12, `pic ip mesuré=${pic}`);
});

test('classique (seed 42) : durée du quench de courant tEnd − tCQ dans [10, 50] ms', () => {
  const { S } = runShot(CLASSIQUE, 42);
  const duree = S.tEnd - S.tCQ;
  assert.ok(duree >= 10 && duree <= 50, `tEnd−tCQ=${duree}`);
});

test('classique (seed 42) : le couplage nourrit le 3/2 — W32 au TQ dans [0.05, 0.20]', () => {
  // Référence mesurée : W32 = 0.1062 au moment du quench thermique (35× la
  // graine). Sans couplage (CPL=0), W32 resterait à W32SEED=0.003 et le TQ
  // n'arriverait que par le seul 2/1 — invariant ajouté car le harnais ne
  // détectait pas la perte du terme de couplage (test par mutation).
  const P = { ...CLASSIQUE };
  const rng = M.mulberry32(42);
  const g = M.makeGauss(rng);
  const S = M.newState(P);
  let w32auTQ = -1;
  while (S.t < M.MP.TMAX && !(S.ended && S.t > S.tEnd + 60)) {
    M.stepModel(S, P, g, DT);
    if (S.tq && w32auTQ < 0) w32auTQ = S.W32;
  }
  assert.ok(w32auTQ >= 0.05 && w32auTQ <= 0.20, `W32 au TQ=${w32auTQ}`);
});

test('sain (seed 42, d0=−0.5) : aucun événement sur 2500 ms', () => {
  const { S, t } = runShot(SAIN, 42);
  assert.equal(S.tLock, -1, `tLock=${S.tLock}`);
  assert.equal(S.tTQ, -1, `tTQ=${S.tTQ}`);
  assert.equal(S.tCQ, -1, `tCQ=${S.tCQ}`);
  assert.equal(S.ended, false);
  assert.ok(t[t.length - 1] >= M.MP.TMAX - DT, 'le tir doit couvrir les 2500 ms');
});

test('sain (seed 42) : T_e (état vrai) toujours > 0.8', () => {
  // Porte sur S.Te, pas sur le canal mesuré : le bruit AR(1) (σ=1.5 %) fait
  // ponctuellement descendre la MESURE sous 0.8 (min observé 0.788) alors que
  // l'état ne descend jamais sous ~0.83.
  const { teEtat } = runShot(SAIN, 42);
  const min = Math.min(...teEtat);
  assert.ok(min > 0.8, `min S.Te=${min}`);
});

test('sain (seed 42) : au moins 80 dents de scie', () => {
  const { dents } = runShot(SAIN, 42);
  assert.ok(dents >= 80, `dents de scie=${dents}`);
});

test('classique (seed 42) : P_rad pique au quench thermique', () => {
  // Références mesurées : pic prad = 1.203 dans [tTQ, tTQ+5 ms] (flash
  // radiatif PRSP), max pré-TQ = 0.237 (montée douce ∝ W21+W32).
  const { t, prad, S } = runShot(CLASSIQUE, 42);
  let picTQ = -1, maxPre = 0;
  for (let i = 0; i < t.length; i++) {
    if (S.tTQ >= 0 && t[i] >= S.tTQ && t[i] <= S.tTQ + 5) picTQ = Math.max(picTQ, prad[i]);
    if (t[i] < S.tTQ && prad[i] > maxPre) maxPre = prad[i];
  }
  assert.ok(picTQ >= 0.8 && picTQ <= 1.6, `pic prad au TQ=${picTQ}`);
  assert.ok(maxPre < 0.5, `prad max pré-TQ=${maxPre}`);
});

test('classique (seed 42) : bouffée de densité au quench thermique', () => {
  // Référence mesurée : ne max = 1.489 (influx NESP au TQ), ne ≈ 1 avant.
  const { ne } = runShot(CLASSIQUE, 42);
  const max = Math.max(...ne);
  assert.ok(max >= 1.2 && max <= 1.8, `ne max=${max}`);
});

test('sain (seed 42) : P_rad bas et densité stable sur tout le tir', () => {
  // Références mesurées : prad max = 0.113, ne ∈ [0.970, 1.026].
  const { prad, ne } = runShot(SAIN, 42);
  assert.ok(Math.max(...prad) < 0.3, `prad max=${Math.max(...prad)}`);
  assert.ok(Math.min(...ne) > 0.9 && Math.max(...ne) < 1.1,
    `ne ∈ [${Math.min(...ne)}, ${Math.max(...ne)}]`);
});

test('reproductibilité : même seed ⇒ séries strictement identiques', () => {
  const a = runShot(CLASSIQUE, 42);
  const b = runShot(CLASSIQUE, 42);
  assert.equal(a.t.length, b.t.length);
  for (let i = 0; i < a.t.length; i++) {
    assert.equal(a.mir[i], b.mir[i], `mir diverge au pas ${i}`);
    assert.equal(a.te[i], b.te[i], `te diverge au pas ${i}`);
    assert.equal(a.ip[i], b.ip[i], `ip diverge au pas ${i}`);
  }
});

test('reproductibilité : seeds différents ⇒ séries différentes', () => {
  const a = runShot(CLASSIQUE, 42);
  const b = runShot(CLASSIQUE, 43);
  const n = Math.min(a.t.length, b.t.length);
  let differe = false;
  for (let i = 0; i < n && !differe; i++) {
    if (a.mir[i] !== b.mir[i] || a.te[i] !== b.te[i] || a.ip[i] !== b.ip[i]) differe = true;
  }
  assert.ok(differe, 'les tirs seed 42 et 43 ne doivent pas être identiques');
});

test('bruit : σ stationnaire du canal te ≈ 0.015 ± 25 % (tir sain, tendance retirée)', () => {
  // Le bruit du canal te est exactement te_mesuré − S.Te (bruit purement
  // observationnel). Fenêtre stationnaire : t ∈ [500, 2500] ms. Théorie AR(1)
  // n' = 0.97·n + sg·g, sg = (1.5/100)·0.2431 ⇒ σ = sg/√(1−0.97²) ≈ 0.0150.
  const { t, te, teEtat } = runShot(SAIN, 42);
  const bruit = [];
  for (let i = 0; i < t.length; i++) {
    if (t[i] >= 500) bruit.push(te[i] - teEtat[i]);
  }
  const moy = bruit.reduce((s, x) => s + x, 0) / bruit.length;
  const variance = bruit.reduce((s, x) => s + (x - moy) * (x - moy), 0) / bruit.length;
  const sigma = Math.sqrt(variance);
  assert.ok(sigma >= 0.015 * 0.75 && sigma <= 0.015 * 1.25, `σ mesuré=${sigma}`);
});
