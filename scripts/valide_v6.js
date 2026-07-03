#!/usr/bin/env node
// Batterie de validation du moteur v6 — à lancer avant toute campagne de
// recherche ou après toute modification physique :
//
//   node scripts/valide_v6.js [--rapide]
//
// Quatre sections, chacune PASS/FAIL, code de sortie non nul si échec :
//   1. Convergence en dt (l'intégrateur d'Euler ne doit pas piloter les
//      résultats) — dynamique isolée des tirages (elm=0).
//   2. Balayage extrême : N tirs AU-DELÀ des plages du générateur, garde
//      anti-NaN et bornes physiques à chaque pas, ordre des événements.
//   3. Propriétés métamorphiques : le bruit de mesure ne change pas la
//      dynamique ; monotonies physiques (d0a, cwf, pheat) ; σ du bruit
//      indépendant du pas.
//   4. Couche diagnostics : finitude et bornes de tous les canaux sur un
//      tir disruptif complet.
'use strict';
const path = require('path');
const M = require(path.join(__dirname, '..', 'model', 'tokamak_model_v6.js'));
const Dg = require(path.join(__dirname, '..', 'model', 'tokamak_diags_v6.js'));

const RAPIDE = process.argv.includes('--rapide');
let nEchecs = 0;

function verdict(nom, ok, detail) {
  console.log((ok ? '  PASS ' : '  FAIL ') + nom + (detail ? ' — ' + detail : ''));
  if (!ok) nEchecs++;
}

function run(P, seed, dtScale, garde) {
  const rng = M.mulberry32(seed), g = M.makeGauss(rng), S = M.newState(P, rng);
  const scale = dtScale || 1;
  while (S.t < M.C.TMAX && !(S.ended && S.t > S.tEnd + 0.05)) {
    M.stepModel(S, P, g, M.chooseDt(S) * scale);
    if (garde && !garde(S)) break;
  }
  return S;
}

// ---------------------------------------------------------------- Section 1
console.log('\n[1] Convergence en dt (dt, dt/2, dt/4 ; elm=0, dynamique déterministe)');
{
  const scenarios = [
    ['classique', { d0a: 1, elm: 0 }],
    ['rapide', { d0a: 2.5, elm: 0 }],
    ['NTM', { d0a: -0.3, bs: 0.5, elm: 0 }],
  ];
  for (const [nom, P] of scenarios) {
    const a = run(P, 42, 1), b = run(P, 42, 0.5), c = run(P, 42, 0.25);
    const derive = (x, y) => Math.abs(x - y) / Math.max(y, 1e-9);
    const dTQ = derive(a.tTQ, c.tTQ), dLock = derive(a.tLock, c.tLock), dEnd = derive(a.tEnd, c.tEnd);
    const cauchy = derive(b.tTQ, c.tTQ) <= derive(a.tTQ, c.tTQ) + 1e-3; // l'erreur doit décroître avec dt
    verdict('convergence ' + nom, dTQ < 0.01 && dLock < 0.01 && dEnd < 0.01 && cauchy,
      `tTQ ${a.tTQ.toFixed(4)} → ${c.tTQ.toFixed(4)} s (dérive ${(100 * dTQ).toFixed(2)} %), ` +
      `tLock dérive ${(100 * dLock).toFixed(2)} %, tEnd dérive ${(100 * dEnd).toFixed(2)} %`);
  }
}

// ---------------------------------------------------------------- Section 2
console.log('\n[2] Balayage extrême (hors plages du générateur), garde anti-NaN et bornes');
{
  const N = RAPIDE ? 100 : 400;
  const meta = M.mulberry32(777);
  const tire = (lo, hi) => lo + meta() * (hi - lo);
  const MJ = M.resolveMachine();
  let nDisrupt = 0, nSainsDisrupt = 0, violations = [], omClampe = 0;
  const wLocks = [], cqNorm = [];
  for (let i = 0; i < N; i++) {
    const sain = i % 2 === 0;
    const P = {
      d0a: sain ? tire(-1.5, -0.2) : tire(0.5, 4),
      f0k: tire(0.5, 8), cwf: tire(0.3, 3),
      pheat: tire(4, 30), n19: tire(1.5, 6), noise: tire(0, 5),
      elm: [0, 0.3, 0.9][i % 3],
      bs: sain ? 0 : [0, 0.3, 0.8][Math.floor(meta() * 3)],
    };
    const seed = 9000 + i;
    const rng = M.mulberry32(seed), g = M.makeGauss(rng), S = M.newState(P, rng);
    let ipPre = -1, tIp80 = -1, tIp20 = -1, omEtaitNul = false;
    while (S.t < M.C.TMAX && !(S.ended && S.t > S.tEnd + 0.05)) {
      M.stepModel(S, P, g, M.chooseDt(S));
      const vals = [S.w, S.w32, S.Om, S.Te0, S.Ip, S.Prad, S.K, S.Ne];
      if (!vals.every(Number.isFinite)) { violations.push({ seed, P, pb: 'non fini', t: S.t }); break; }
      if (S.Te0 < -1e-6 || S.Te0 > 40) { violations.push({ seed, P, pb: 'Te0 hors bornes', Te0: S.Te0 }); break; }
      if (S.w > M.C.WSAT * MJ.A * 1.0001 || S.w < 0) { violations.push({ seed, P, pb: 'w hors bornes', w: S.w }); break; }
      if (S.Ip < -1 || S.Ip > 1.2 * MJ.IP0) { violations.push({ seed, P, pb: 'Ip hors bornes', Ip: S.Ip }); break; }
      if (S.Om < 0 || S.Om > 1.05 * S.Om0) { violations.push({ seed, P, pb: 'Om hors bornes', Om: S.Om }); break; }
      if (S.Om === 0 && !S.locked) omEtaitNul = true;
      if (S.cq && ipPre < 0) ipPre = S.Ip;
      if (ipPre > 0) {
        if (tIp80 < 0 && S.Ip <= 0.8 * ipPre) tIp80 = S.t;
        if (tIp20 < 0 && S.Ip <= 0.2 * ipPre) tIp20 = S.t;
      }
    }
    if (omEtaitNul) omClampe++;
    const evOk = (S.tCQ < 0 || (S.tTQ >= 0 && S.tTQ < S.tCQ)) &&
                 (S.tEnd < 0 || (S.tCQ >= 0 && S.tCQ < S.tEnd));
    if (!evOk) violations.push({ seed, P, pb: 'ordre des événements' });
    if (S.tTQ >= 0) nDisrupt++;
    if (sain && S.tTQ >= 0) nSainsDisrupt++;
    if (S.wAtLock > 0) wLocks.push(S.wAtLock * 100);
    if (tIp80 > 0 && tIp20 > 0) cqNorm.push((tIp20 - tIp80) / (Math.PI * MJ.A * MJ.A * MJ.KAPPA) * 1e3);
  }
  verdict('aucune violation de bornes / NaN sur ' + N + ' tirs', violations.length === 0,
    violations.length ? JSON.stringify(violations.slice(0, 3)) : 'RAS');
  verdict('aucun tir sain (d0a<0, bs=0) ne disrompt', nSainsDisrupt === 0, nSainsDisrupt + ' cas');
  wLocks.sort((a, b) => a - b); cqNorm.sort((a, b) => a - b);
  const q = (arr, p) => arr[Math.floor(p * (arr.length - 1))];
  verdict('largeurs au verrouillage dans [1, 8] cm (médiane ' + q(wLocks, 0.5).toFixed(1) + ')',
    wLocks.length > 0 && q(wLocks, 0) >= 1 && q(wLocks, 1) <= 8,
    `min ${q(wLocks, 0).toFixed(1)} / max ${q(wLocks, 1).toFixed(1)} cm sur ${wLocks.length} verrouillages`);
  verdict('quench 80→20 % normalisé dans [1.7, 30] ms/m²',
    cqNorm.length > 0 && q(cqNorm, 0) >= 1.7 && q(cqNorm, 1) <= 30,
    `min ${q(cqNorm, 0).toFixed(1)} / méd ${q(cqNorm, 0.5).toFixed(1)} / max ${q(cqNorm, 1).toFixed(1)} ms/m² sur ${cqNorm.length} quenchs`);
  console.log(`  info : ${nDisrupt}/${N} disruptions, clamp Ω≥0 sollicité sur ${omClampe} tir(s)`);
}

// ---------------------------------------------------------------- Section 3
console.log('\n[3] Propriétés métamorphiques');
{
  // (a) le bruit de mesure ne change pas la dynamique (mêmes tirages, autre σ)
  const a = run({ d0a: 1, noise: 0 }, 42, 1), b = run({ d0a: 1, noise: 5 }, 42, 1);
  verdict('bruit sans effet sur la dynamique (événements identiques à σ=0 et σ=5 %)',
    a.tLock === b.tLock && a.tTQ === b.tTQ && a.tCQ === b.tCQ && a.tEnd === b.tEnd,
    `tTQ ${a.tTQ} vs ${b.tTQ}`);
  // (b) plus d'instabilité ⇒ disruption plus tôt
  const grille = [0.6, 1.0, 1.5, 2.2, 3.0].map(d0a => run({ d0a, elm: 0 }, 42, 1).tTQ);
  verdict('monotonie : tTQ décroît avec d0a', grille.every((t, i) => i === 0 || t < grille[i - 1]),
    grille.map(t => t.toFixed(2)).join(' > '));
  // (c) plus de freinage ⇒ verrouillage plus tôt et îlot plus petit
  const cw = [0.5, 1, 2].map(cwf => run({ d0a: 1, cwf, elm: 0 }, 42, 1));
  verdict('monotonie : tLock et w_lock décroissent avec le freinage cwf',
    cw[0].tLock > cw[1].tLock && cw[1].tLock > cw[2].tLock &&
    cw[0].wAtLock > cw[1].wAtLock && cw[1].wAtLock > cw[2].wAtLock,
    cw.map(s => `${s.tLock.toFixed(2)}s/${(s.wAtLock * 100).toFixed(1)}cm`).join(' , '));
  // (d) plus de chauffage ⇒ plasma plus chaud (tir sain)
  const teMax = [6, 10, 16].map(pheat => {
    const rng = M.mulberry32(42), g = M.makeGauss(rng), S = M.newState({ d0a: -0.5, pheat, elm: 0 }, rng);
    let mx = 0;
    while (S.t < 3) { M.stepModel(S, { d0a: -0.5, pheat, elm: 0 }, g, M.chooseDt(S)); if (S.Te0 > mx) mx = S.Te0; }
    return mx;
  });
  verdict('monotonie : Te max croît avec P_chauffage', teMax[0] < teMax[1] && teMax[1] < teMax[2],
    teMax.map(t => t.toFixed(1) + ' keV').join(' < '));
  // (e) σ stationnaire du bruit indépendant du pas d'intégration (AR(1) en dt)
  const sigma = (scale) => {
    const P = { d0a: -0.5, elm: 0 };
    const rng = M.mulberry32(42), g = M.makeGauss(rng), S = M.newState(P, rng);
    const v = [];
    while (S.t < 2) { const m = M.stepModel(S, P, g, 1e-4 * scale); if (S.t > 0.5) v.push(m.te - S.Te0); }
    const mu = v.reduce((s, x) => s + x, 0) / v.length;
    return Math.sqrt(v.reduce((s, x) => s + (x - mu) ** 2, 0) / v.length);
  };
  const s1 = sigma(1), s2 = sigma(0.5);
  verdict('σ du bruit indépendant du pas (dt vs dt/2)', Math.abs(s1 - s2) / s1 < 0.1,
    `${s1.toFixed(4)} vs ${s2.toFixed(4)} keV`);
}

// ---------------------------------------------------------------- Section 4
console.log('\n[4] Couche diagnostics : finitude et bornes sur un tir disruptif complet');
{
  const P = { d0a: 1.5, f0k: 2 };
  const rng = M.mulberry32(4242), g = M.makeGauss(rng), S = M.newState(P, rng);
  const p = M.defP(P), D = Dg.newDiag(P, S);
  let tS = 0, nRows = 0, pb = null, maxCoil = 0, maxSad = 0, maxBolo = 0;
  while (S.t < M.C.TMAX && !(S.ended && S.t > S.tEnd + 0.05)) {
    M.stepModel(S, P, g, M.chooseDt(S));
    if (S.t >= tS) {
      const row = Dg.sampleDiag(D, S, p, g, 1e-4);
      nRows++;
      if (!row.every(Number.isFinite)) { pb = 'valeur non finie à t=' + S.t; break; }
      for (let k = 1; k <= 8; k++) maxCoil = Math.max(maxCoil, Math.abs(row[k]));
      for (let k = 9; k <= 12; k++) maxSad = Math.max(maxSad, Math.abs(row[k]));
      maxBolo = Math.max(maxBolo, row[19]);
      tS += 1e-4;
    }
  }
  verdict('tous les canaux finis sur ' + nRows + ' échantillons', pb === null, pb || 'RAS');
  verdict('bornes physiques des capteurs', maxCoil < 500 && maxSad < 0.1 && maxBolo < 20000,
    `|dB/dt| max ${maxCoil.toFixed(1)} T/s, |Br| max ${(maxSad * 1e3).toFixed(2)} mT, bolo max ${maxBolo.toFixed(0)} MW`);
}

console.log('\n' + (nEchecs === 0 ? 'VALIDATION COMPLÈTE : tout est vert.'
                                  : 'ÉCHECS : ' + nEchecs + ' vérification(s) en défaut.'));
process.exitCode = nEchecs === 0 ? 0 : 1;
