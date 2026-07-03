// Tests du modèle v6 — ancrés dans la littérature, pas dans le modèle.
// Chaque plage encode un chiffre ou une bande PUBLIÉE (voir docs/V6_PHYSIQUE.md) :
// le test échoue si la physique dimensionnée dérive de l'expérience, pas si un
// détail numérique interne bouge. Références : ITER Physics Basis 1999 (IPB98,
// plancher CQ 1.7 ms/m²), de Vries NF 51 (2011) (précurseurs JET), Sweeney
// NF 57 (2017) (modes verrouillés), Wesson « Tokamaks » (Spitzer, L_p, TQ).
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const M = require('../model/tokamak_model_v6.js');

const SAIN = { d0a: -0.5 };
const CLASSIQUE = { d0a: 1 };
const MJ = M.resolveMachine('jet-like');

// Un tir complet à pas adaptatif (chooseDt), avec instrumentation.
function runShot(P, seed) {
  const rng = M.mulberry32(seed), g = M.makeGauss(rng), S = M.newState(P, rng);
  const o = {
    t: [], mir: [], te: [], ip: [], prad: [], ne: [], teEtat: [],
    dents: 0, elms: 0, maxMirPh1: 0, maxPradMW: 0, pradPreTQ: 0,
    ipPeakMA: 0, teMinApres05: 99, teMaxApres05: 0,
    tqT80: -1, tqT20: -1, cqT80: -1, cqT20: -1,
  };
  o.omAtLock = -1;
  let tePrec = S.Te0, elmPrec = 0, tePreTQ = -1, ipPreCQ = -1, tSample = 0;
  while (S.t < M.C.TMAX && !(S.ended && S.t > S.tEnd + 0.1)) {
    const dt = M.chooseDt(S);
    const meas = M.stepModel(S, P, g, dt);
    if (S.locked && o.omAtLock < 0) o.omAtLock = S.Om / S.Om0;
    if (S.t >= tSample) { // échantillonnage 1 kHz
      o.t.push(S.t); o.mir.push(meas.mir); o.te.push(meas.te); o.ip.push(meas.ip);
      o.prad.push(meas.prad); o.ne.push(meas.ne); o.teEtat.push(S.Te0);
      tSample += 1e-3;
    }
    if (!S.tq) {
      if (S.Te0 < tePrec * 0.9) o.dents++;
      if (S.t > 0.5) {
        if (S.Te0 < o.teMinApres05) o.teMinApres05 = S.Te0;
        if (S.Te0 > o.teMaxApres05) o.teMaxApres05 = S.Te0;
      }
      if (meas.prad > o.pradPreTQ) o.pradPreTQ = meas.prad;
    }
    if (S.elmB > elmPrec + 0.3) o.elms++;
    elmPrec = S.elmB;
    if (S.phaseId === 1 && Math.abs(meas.mir) > o.maxMirPh1) o.maxMirPh1 = Math.abs(meas.mir);
    if (meas.prad > o.maxPradMW) o.maxPradMW = meas.prad;
    if (S.Ip / 1e6 > o.ipPeakMA) o.ipPeakMA = S.Ip / 1e6;
    // Chute de Te 80 % → 20 % au quench thermique
    if (S.tq && tePreTQ < 0) tePreTQ = tePrec;
    if (tePreTQ > 0) {
      if (o.tqT80 < 0 && S.Te0 <= 0.8 * tePreTQ) o.tqT80 = S.t;
      if (o.tqT20 < 0 && S.Te0 <= 0.2 * tePreTQ) o.tqT20 = S.t;
    }
    // Chute d'Ip 80 % → 20 % au quench de courant
    if (S.cq && ipPreCQ < 0) ipPreCQ = S.Ip;
    if (ipPreCQ > 0) {
      if (o.cqT80 < 0 && S.Ip <= 0.8 * ipPreCQ) o.cqT80 = S.t;
      if (o.cqT20 < 0 && S.Ip <= 0.2 * ipPreCQ) o.cqT20 = S.t;
    }
    tePrec = S.Te0;
  }
  o.S = S;
  return o;
}

// Les tirs de référence sont coûteux (~100 k pas) : un seul run par scénario.
const cache = new Map();
function shot(P, seed) {
  const k = JSON.stringify([P, seed]);
  if (!cache.has(k)) cache.set(k, runShot(P, seed));
  return cache.get(k);
}

test('Spitzer : η(2 keV) ≈ 1e-8 Ω·m (Wesson)', () => {
  const eta = M.etaSpitzer(2, 1.7, 1);
  assert.ok(eta >= 8e-9 && eta <= 1.2e-8, `η=${eta}`);
});

test('IPB98(y,2) : τ_E JET-like (2.5 MA, n̄=3, 10 MW) dans [0.2, 0.35] s', () => {
  const tau = M.tauE98(MJ, 2.5e6, 3, 1e7);
  assert.ok(tau >= 0.2 && tau <= 0.35, `τ_E=${tau}`);
  // Dépendances du scaling : favorable en courant, dégradé en puissance.
  assert.ok(M.tauE98(MJ, 3.5e6, 3, 1e7) > tau, 'τ_E doit croître avec Ip');
  assert.ok(M.tauE98(MJ, 2.5e6, 3, 2e7) < tau, 'τ_E doit décroître avec P');
});

test('IPB98(y,2) : le preset iter-like retrouve la prédiction publiée τ_E ≈ 3.7 s', () => {
  // Scénario Q=10 : Ip = 15 MA, n̄ ≈ 10×10¹⁹ m⁻³, P_pertes ≈ 87 MW
  // (ITER Physics Basis : τ_E prédit ≈ 3.7 s).
  const MI = M.resolveMachine('iter-like');
  const tau = M.tauE98(MI, 15e6, 10, 87e6);
  assert.ok(tau >= 3.0 && tau <= 4.5, `τ_E ITER=${tau}`);
});

test('machine : presets et surcharges via un objet de config unique', () => {
  // Un nom inconnu doit échouer clairement.
  assert.throws(() => M.resolveMachine('sparc'), /machine inconnue/);
  // Une surcharge partielle recalcule les grandeurs dérivées.
  const petit = M.resolveMachine({ A: 0.5 });
  assert.equal(petit.R0, MJ.R0);
  assert.ok(Math.abs(petit.VOL - MJ.VOL / 4) < 1e-9, 'VOL ∝ a²');
  assert.ok(Math.abs(petit.DRS - MJ.DRS / 2) < 1e-12, 'Δr_s ∝ a');
  // La machine du tir est bien celle demandée.
  const S = M.newState({ machine: 'iter-like' });
  assert.equal(S.Ip, 15e6);
});

test('sain (d0a=−0.5, seed 42) : 10 s sans événement, Te0 JET-like [4, 8] keV', () => {
  const o = shot(SAIN, 42);
  assert.equal(o.S.tLock, -1);
  assert.equal(o.S.tTQ, -1);
  assert.equal(o.S.tCQ, -1);
  assert.ok(o.teMinApres05 >= 4 && o.teMaxApres05 <= 8,
    `Te0 ∈ [${o.teMinApres05}, ${o.teMaxApres05}] keV`);
});

test('sain : dents de scie et ELM à des cadences réalistes', () => {
  const o = shot(SAIN, 42);
  // Période de dents ~0.12 s ±10 % ⇒ 75-93 crashs sur 10 s (JET : 0.05-0.3 s).
  assert.ok(o.dents >= 60 && o.dents <= 110, `dents=${o.dents}`);
  // ELM type-I ~20-50 Hz.
  const freq = o.elms / M.C.TMAX;
  assert.ok(freq >= 20 && freq <= 50, `f_ELM=${freq} Hz`);
});

test('classique (d0a=1, seed 42) : chaîne complète ordonnée onset → lock → TQ → CQ → fin', () => {
  const { S } = shot(CLASSIQUE, 42);
  assert.ok(S.tOnset >= 0 && S.tOnset < S.tLock, `tOnset=${S.tOnset}`);
  assert.ok(S.tLock < S.tTQ && S.tTQ < S.tCQ && S.tCQ < S.tEnd,
    `${S.tLock} < ${S.tTQ} < ${S.tCQ} < ${S.tEnd}`);
});

test('classique : verrouillage à une largeur d’îlot de quelques cm (Sweeney, 3 seeds)', () => {
  for (const seed of [42, 43, 44]) {
    const { S } = shot(CLASSIQUE, seed);
    assert.ok(S.wAtLock >= 0.015 && S.wAtLock <= 0.06, `seed ${seed} : w au lock=${S.wAtLock} m`);
  }
});

test('verrouillage : le drapeau locked marque un vrai effondrement de la rotation', () => {
  // Vérité terrain des datasets ML : au moment où locked passe à vrai, la
  // rotation doit s'être réellement effondrée (Ω ≤ 10 % de Ω0), pas seulement
  // ralentie. Ajouté après test par mutation : LOCKFRAC 0.05→0.3 passait
  // inaperçu car la branche verrouillée écrase Ω ensuite de toute façon.
  for (const seed of [42, 43, 44]) {
    const o = shot(CLASSIQUE, seed);
    assert.ok(o.omAtLock >= 0 && o.omAtLock <= 0.1,
      `seed ${seed} : Ω/Ω0 au verrouillage=${o.omAtLock}`);
  }
});

test('calibration du couple de paroi : w_lock de référence dans [2.5, 4.2] cm', () => {
  // Plus étroit que la bande littérature [1.5, 6] cm du test Sweeney : ce
  // test ÉPINGLE la calibration de C.CW sur le scénario de référence (un
  // CW divisé par 4 donne w_lock ≈ 4.7 cm, indétectable par la bande large —
  // test par mutation). Re-calibrer CW exige de mettre à jour ces bornes,
  // avec justification dans le commit.
  for (const seed of [42, 43, 44]) {
    const { S } = shot(CLASSIQUE, seed);
    assert.ok(S.wAtLock >= 0.025 && S.wAtLock <= 0.042,
      `seed ${seed} : w au lock=${S.wAtLock} m`);
  }
});

test('classique : durées de précurseur à l’échelle JET (de Vries)', () => {
  const { S } = shot(CLASSIQUE, 42);
  const verrouille = S.tTQ - S.tLock, total = S.tTQ - S.tOnset;
  assert.ok(verrouille >= 0.05 && verrouille <= 3, `phase verrouillée=${verrouille} s`);
  assert.ok(total >= 0.3 && total <= 8, `précurseur total=${total} s`);
});

test('classique : quench thermique en 0.05–1.5 ms (chute Te 80→20 %, 3 seeds)', () => {
  for (const seed of [42, 43, 44]) {
    const o = shot(CLASSIQUE, seed);
    const duree = (o.tqT20 - o.tqT80) * 1e3;
    assert.ok(o.tqT80 > 0 && duree >= 0.05 && duree <= 1.5, `seed ${seed} : TQ 80→20=${duree} ms`);
  }
});

test('classique : pic d’Ip de +2 à +10 % (aplatissement du profil)', () => {
  const o = shot(CLASSIQUE, 42);
  const pic = o.ipPeakMA / (MJ.IP0 / 1e6);
  assert.ok(pic >= 1.02 && pic <= 1.10, `pic Ip=${pic}×Ip0`);
});

test('classique : quench de courant 80→20 % normalisé dans [1.7, 30] ms/m² (ITER/JET, 3 seeds)', () => {
  for (const seed of [42, 43, 44]) {
    const o = shot(CLASSIQUE, seed);
    const norme = (o.cqT20 - o.cqT80) / (Math.PI * MJ.A * MJ.A * MJ.KAPPA) * 1e3;
    assert.ok(o.cqT80 > 0 && norme >= 1.7 && norme <= 30, `seed ${seed} : CQ normalisé=${norme} ms/m²`);
    const totale = o.S.tEnd - o.S.tCQ;
    assert.ok(totale >= 0.02 && totale <= 0.15, `seed ${seed} : durée CQ totale=${totale} s`);
  }
});

test('fiabilité : convergence en dt — les événements ne dépendent pas du pas', () => {
  // Dynamique isolée des tirages (elm=0) ; dt et dt/2 doivent donner les
  // mêmes temps d'événements à mieux que 1 % (validation complète : dt/4
  // dans scripts/valide_v6.js).
  const P = { d0a: 1, elm: 0 };
  const run = (scale) => {
    const rng = M.mulberry32(42), g = M.makeGauss(rng), S = M.newState(P, rng);
    while (S.t < M.C.TMAX && !(S.ended && S.t > S.tEnd + 0.05)) {
      M.stepModel(S, P, g, M.chooseDt(S) * scale);
    }
    return S;
  };
  const a = run(1), b = run(0.5);
  for (const ev of ['tLock', 'tTQ', 'tEnd']) {
    assert.ok(Math.abs(a[ev] - b[ev]) / b[ev] < 0.01, `${ev} : ${a[ev]} vs ${b[ev]}`);
  }
});

test('fiabilité : le bruit de mesure ne modifie jamais la dynamique', () => {
  // Le bruit est purement observationnel : mêmes tirages, seule l'échelle σ
  // change ⇒ les événements doivent être STRICTEMENT identiques.
  const a = shot({ d0a: 1, noise: 0 }, 42).S;
  const b = shot({ d0a: 1, noise: 5 }, 42).S;
  assert.equal(a.tLock, b.tLock);
  assert.equal(a.tTQ, b.tTQ);
  assert.equal(a.tCQ, b.tCQ);
  assert.equal(a.tEnd, b.tEnd);
});

test('fiabilité : tTQ décroît strictement avec l’instabilité d0a', () => {
  const ts = [0.8, 1.4, 2.4].map(d0a => shot({ d0a, elm: 0 }, 42).S.tTQ);
  assert.ok(ts[0] > ts[1] && ts[1] > ts[2], ts.join(' > '));
});

test('CQ : la durée varie d’un tir à l’autre (plateau post-quench 6–16 eV)', () => {
  // Les distributions expérimentales de temps de quench s'étalent sur
  // presque un ordre de grandeur ; un modèle à durée unique serait irréaliste.
  const a = shot({ d0a: 1.5 }, 100);
  const b = shot({ d0a: 1.5 }, 101);
  const da = a.S.tEnd - a.S.tCQ, db = b.S.tEnd - b.S.tCQ;
  assert.ok(da > 0 && db > 0);
  assert.ok(Math.abs(da - db) / Math.max(da, db) > 0.03,
    `durées CQ trop semblables : ${da} vs ${db} s`);
});

test('classique : précurseur Mirnov de l’ordre du T/s (bobines réelles)', () => {
  const o = shot(CLASSIQUE, 42);
  assert.ok(o.maxMirPh1 >= 0.5 && o.maxMirPh1 <= 50, `max |dB/dt| phase 1=${o.maxMirPh1} T/s`);
});

test('classique : flash radiatif du TQ à l’échelle GW, rayonnement calme sinon', () => {
  const o = shot(CLASSIQUE, 42);
  assert.ok(o.maxPradMW >= 500 && o.maxPradMW <= 20000, `flash=${o.maxPradMW} MW`);
  assert.ok(o.pradPreTQ >= 1 && o.pradPreTQ <= 25, `P_rad pré-TQ=${o.pradPreTQ} MW`);
});

test('NTM : métastabilité — stable classiquement, disruptif avec bootstrap', () => {
  const sans = shot({ d0a: -0.3, bs: 0 }, 42);
  const avec = shot({ d0a: -0.3, bs: 0.5 }, 42);
  assert.equal(sans.S.tTQ, -1, 'sans bootstrap : pas de disruption');
  assert.ok(avec.S.tTQ > 0, 'avec bootstrap : l’îlot doit croître et disrompre');
});

test('reproductibilité : même seed ⇒ identique, seeds ≠ ⇒ différents', () => {
  const a = runShot(CLASSIQUE, 42);
  const b = runShot(CLASSIQUE, 42);
  assert.deepEqual(a.te, b.te);
  assert.deepEqual(a.mir, b.mir);
  const c = runShot(CLASSIQUE, 43);
  const n = Math.min(a.te.length, c.te.length);
  let differe = false;
  for (let i = 0; i < n && !differe; i++) if (a.te[i] !== c.te[i]) differe = true;
  assert.ok(differe);
});

test('bruit : σ stationnaire du canal te ≈ noise%×5 keV ± 25 % (AR(1) correct en dt)', () => {
  const o = shot(SAIN, 42);
  const bruit = [];
  for (let i = 0; i < o.t.length; i++) {
    if (o.t[i] >= 2) bruit.push(o.te[i] - o.teEtat[i]);
  }
  const moy = bruit.reduce((s, x) => s + x, 0) / bruit.length;
  const sigma = Math.sqrt(bruit.reduce((s, x) => s + (x - moy) ** 2, 0) / bruit.length);
  const attendu = 0.015 * 5;
  assert.ok(sigma >= attendu * 0.75 && sigma <= attendu * 1.25, `σ=${sigma} keV`);
});

test('chooseDt : pas fin (10 µs) pendant les quenches, grossier (100 µs) sinon', () => {
  const rng = M.mulberry32(1), S = M.newState({}, rng);
  assert.equal(M.chooseDt(S), 1e-4);
  S.tq = true;
  assert.equal(M.chooseDt(S), 1e-5);
  S.ended = true;
  assert.equal(M.chooseDt(S), 1e-4);
});
