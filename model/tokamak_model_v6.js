/* Modèle tokamak v6 — 0D en unités SI, machine JET-like, calibré littérature.
   Voir docs/V6_PHYSIQUE.md pour les équations et les références. Ce moteur est
   distinct de la v5 pédagogique (model/tokamak_model.js), qui reste intacte.
   Unités : t en s, w en m, Ω en rad/s, Te0 en keV, Ip en A, Prad en W.
   Double export : window.TokamakModelV6 en navigateur, module.exports en Node.
   Toute la stochasticité passe par le rng/gauss fournis (reproductible). */
(function (root, factory) {
  if (typeof module === 'object' && typeof module.exports === 'object') {
    module.exports = factory();
  } else {
    root.TokamakModelV6 = factory();
  }
})(typeof self !== 'undefined' ? self : this, function () {
/*MODEL-V6-BEGIN*/
var C = {
  MU0: 4e-7 * Math.PI,
  // --- Machine JET-like (section circulaire)
  R0: 3.0, A: 1.0, B0: 3.0, KAPPA: 1.0, MION: 2, ZEFF: 1.7, LNL: 17,
  IP0: 2.5e6,                    // A
  RS21: 0.75, RS32: 0.60,        // rayons q=2 et q=1.5 (fractions de a)
  // --- Profils (figés en 0D)
  TERSFRAC: 0.4,                 // Te(r_s)/Te0
  PEAK: 2.5,                     // Te0 / <T>
  NEO: 2.5,                      // correction néoclassique de la résistivité
  LI: 1.0,                       // inductance interne
  // --- Îlots (fractions de a sauf mention)
  WSEED: 0.005, WSAT: 0.25, W32SAT: 0.20,
  D32A: -0.5, CPL: 8, LBOOST: 1.5,
  WDBS: 0.02,                    // largeur de coupure du terme bootstrap
  // --- Rotation / paroi
  TAUW: 0.005,                   // temps résistif de la chambre (s)
  CW: 4.0e6,                     // couple de paroi (calibré : lock à w ~ 2-5 cm)
  LOCKFRAC: 0.05, TAULOCK: 0.010,
  // --- Énergie (IPB98(y,2), H98 = 1)
  H98: 1.0, FDEGK: 1.5, FDEGMIN: 0.15,
  // --- Dents de scie / ELM (s)
  SAWP: 0.12, SAWJIT: 0.1, SAWDROP: 0.85, SAWTE: 2.0,
  ELMP: 0.030, ELMJIT: 0.25, ELMTAU: 0.002, ELMDROP: 0.02, ELMTE: 3.0,
  ELMFHZ: 3000, ELMMIR: 1.5, ELMPRAD: 3e6,
  // --- Quenches
  KCRIT: 1.0, TQTAU: 3e-4, TQTE: 0.05, CQTE: 0.10, TECQ: 0.010,
  ZEFFCQ: 3, IPSPIKE: 0.05, ENDIP: 0.05,
  // --- Mesures
  CMIR: 0.05, MIR32: 0.35,
  PRADBASE: 0.3,                 // fraction rayonnée de P_chauffage hors quench
  PRADW: 30e6,                   // W de rayonnement par (w+w32)/a (bord froid)
  NECQ: 0.4, NEINFLUX: 0.5, NETAU: 0.25, NESPT: 0.006,
  NOISETAU: 1.5e-3,              // temps de corrélation du bruit AR(1) (s)
  TMAX: 10.0,
};
// Grandeurs dérivées de la machine
C.VOL = 2 * Math.PI * Math.PI * C.R0 * C.A * C.A * C.KAPPA;          // m^3
C.KW = 3 * 1.602e-16 * 1e19 * C.VOL / C.PEAK;                        // J / (keV·n19)
C.LP = C.MU0 * C.R0 * (Math.log(8 * C.R0 / C.A) - 2 + C.LI / 2);     // H
C.DRS = (C.RS21 - C.RS32) * C.A;                                     // m
C.EPSS = C.RS21 * C.A / C.R0;                                        // ε à r_s

function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;var t=Math.imul(a^a>>>15,1|a);
  t=t+Math.imul(t^t>>>7,61|t)^t;return ((t^t>>>14)>>>0)/4294967296;};}
function makeGauss(rng){var spare=null;return function(){
  if(spare!==null){var s=spare;spare=null;return s;}
  var u=0,v=0;while(u===0)u=rng();while(v===0)v=rng();
  var m=Math.sqrt(-2*Math.log(u)),th=2*Math.PI*v;
  spare=m*Math.sin(th);return m*Math.cos(th);};}

// Résistivité de Spitzer parallèle avec correction néoclassique (Ω·m).
function etaSpitzer(TeKeV, zeff, neo) {
  return (neo || 1) * 1.65e-9 * C.LNL * (zeff / 1.7) /
         Math.pow(Math.max(TeKeV, 0.005), 1.5);
}

// Scaling de confinement IPB98(y,2) — Ip en A, P en W, n en 1e19 m^-3.
function tauE98(IpA, n19, PW) {
  var ipMA = Math.max(IpA, 1e4) / 1e6, pMW = Math.max(PW, 1e5) / 1e6;
  return C.H98 * 0.0562 *
    Math.pow(ipMA, 0.93) * Math.pow(C.B0, 0.15) * Math.pow(n19, 0.41) *
    Math.pow(pMW, -0.69) * Math.pow(C.R0, 1.97) *
    Math.pow(C.A / C.R0, 0.58) * Math.pow(C.KAPPA, 0.78) * Math.pow(C.MION, 0.19);
}

// Beta poloïdal (profil figé) pour le terme bootstrap.
function betaP(Te0, n19, IpA) {
  var pMoy = (2 / 3) * C.KW * n19 * Te0 / C.VOL;          // Pa
  var bTheta = C.MU0 * Math.max(IpA, 1e4) / (2 * Math.PI * C.A);
  return 2 * C.MU0 * pMoy / (bTheta * bTheta);
}

// P : { d0a, f0k, cwf, bs, pheat, n19, noise, elm } — défauts JET-like.
function defP(P) {
  return {
    d0a: P.d0a === undefined ? 1.0 : P.d0a,      // Δ'0·a (sans dimension)
    f0k: P.f0k === undefined ? 5.0 : P.f0k,      // rotation initiale (kHz)
    cwf: P.cwf === undefined ? 1.0 : P.cwf,      // facteur de freinage paroi
    bs: P.bs === undefined ? 0 : P.bs,           // amplitude bootstrap (NTM), 0 = off
    pheat: P.pheat === undefined ? 10 : P.pheat, // MW
    n19: P.n19 === undefined ? 3 : P.n19,        // 1e19 m^-3
    noise: P.noise === undefined ? 1.5 : P.noise,// % (échelle par canal)
    elm: P.elm === undefined ? 0.6 : P.elm,      // amplitude ELM, 0 = off
  };
}

function newState(P, rng) {
  var p = defP(P);
  var om0 = 2 * Math.PI * p.f0k * 1e3;
  return { t: 0,
    w: C.WSEED * C.A, w32: C.WSEED * C.A * 0.6,
    Om: om0, Om0: om0,
    Te0: 3.0,                                    // démarre froid, monte vers l'équilibre
    Ip: C.IP0, TeCQ: -1,
    locked: false, tq: false, cq: false, ended: false,
    sawT: 0, sawP: C.SAWP * (rng ? (1 - C.SAWJIT + 2 * C.SAWJIT * rng()) : 1),
    elmB: 0, elmT: C.ELMP, elmPh: 0,
    rotPh: 0, K: 0, phaseId: 0,
    Prad: 0, Ne: 0, prFlash: 0, neSpike: 0,
    nM: 0, nT: 0, nI: 0, nP: 0, nN: 0,
    tLock: -1, tTQ: -1, tCQ: -1, tEnd: -1, tOnset: -1, wAtLock: -1 };
}

// Bruit AR(1) correct en dt : ρ = exp(−dt/τn), variance stationnaire = σ².
function stepNoise(S, p, g, dt) {
  var rho = Math.exp(-dt / C.NOISETAU);
  var q = Math.sqrt(1 - rho * rho) * (p.noise / 100);
  S.nM = rho * S.nM + q * g();   // × échelle mirnov (T/s) au retour
  S.nT = rho * S.nT + q * g();   // × échelle Te (keV)
  S.nI = rho * S.nI + q * g();   // × échelle Ip (MA)
  S.nP = rho * S.nP + q * g();   // × échelle Prad (MW)
  S.nN = rho * S.nN + q * g();   // × échelle ne (1e19)
}

function measures(S, p) {
  var mir = 0;
  if (!S.ended) {
    var wa = S.w / C.A, w32a = S.w32 / C.A;
    mir = C.CMIR * C.B0 * wa * wa * S.Om * Math.sin(S.rotPh)
        + C.MIR32 * C.CMIR * C.B0 * w32a * w32a * S.Om * Math.sin(1.5 * S.rotPh + 1.1)
        + C.ELMMIR * S.elmB * Math.sin(S.elmPh);
  }
  return {
    mir: mir + S.nM * 1.0,                        // T/s (bruit : échelle 1 T/s)
    te: S.Te0 + S.nT * 5.0,                       // keV (échelle 5 keV)
    ip: S.Ip / 1e6 + S.nI * 2.5,                  // MA (échelle Ip0)
    prad: (S.Prad + S.prFlash) / 1e6 + S.nP * 3,  // MW (échelle 3 MW)
    ne: p.n19 * (S.Ne + 1) + S.neSpike + S.nN * 0.5, // 1e19 m^-3
  };
}

function stepModel(S, P, g, dt) {
  var p = defP(P);
  if (S.ended) {
    S.t += dt;
    S.Te0 += (0 - S.Te0) * dt / 0.05;
    S.Ip += (0 - S.Ip) * dt / 0.05;
    S.Prad += (0 - S.Prad) * dt / 0.05;
    S.prFlash *= Math.exp(-dt / C.TQTAU);
    S.neSpike *= Math.exp(-dt / C.NESPT);
    S.Ne += (-(1 - C.NECQ) - S.Ne) * dt / C.NETAU;
    stepNoise(S, p, g, dt);
    return measures(S, p);
  }

  // --- Rutherford (classique + verrouillage + bootstrap optionnel + couplage)
  var teRs = C.TERSFRAC * S.Te0;
  var gamma = 1.22 * etaSpitzer(teRs, C.ZEFF, C.NEO) / C.MU0;   // m²/s
  var wa = S.w / C.A;
  var dEff = (p.d0a * (1 - S.w / (C.WSAT * C.A)) + (S.locked ? C.LBOOST : 0)) / C.A;
  if (p.bs > 0) {
    var bp = betaP(S.Te0, p.n19, S.Ip);
    dEff += p.bs * Math.sqrt(C.EPSS) * bp * wa / (wa * wa + C.WDBS * C.WDBS) / C.A;
  }
  S.w += gamma * dEff * dt;
  if (S.w < C.WSEED * C.A) S.w = C.WSEED * C.A;
  if (S.w > C.WSAT * C.A) S.w = C.WSAT * C.A;
  if (S.tOnset < 0 && S.w >= 0.01 * C.A) S.tOnset = S.t;   // îlot détectable (1 cm)
  var d32 = (C.D32A + C.CPL * wa) * (1 - S.w32 / (C.W32SAT * C.A)) / C.A;
  S.w32 += gamma * d32 * dt;
  if (S.w32 < C.WSEED * C.A * 0.6) S.w32 = C.WSEED * C.A * 0.6;
  if (S.w32 > C.W32SAT * C.A) S.w32 = C.W32SAT * C.A;

  // --- Rotation : restauration visqueuse vs couple de paroi (∝ w⁴, emballement)
  var tauE = tauE98(S.Ip, p.n19, p.pheat * 1e6);
  if (!S.locked) {
    var brake = (C.CW * p.cwf / C.TAUW) * Math.pow(wa, 4) * S.Om / (1 + S.Om * C.TAUW);
    S.Om += ((S.Om0 - S.Om) / tauE - brake) * dt;
    if (S.Om < C.LOCKFRAC * S.Om0) { S.locked = true; S.tLock = S.t; S.wAtLock = S.w; }
  } else {
    S.Om += (0 - S.Om) * dt / C.TAULOCK;
  }

  // --- Recouvrement (Chirikov) → quench thermique
  S.K = (S.w + S.w32) / C.DRS;
  if (!S.tq && S.K >= C.KCRIT) {
    S.tq = true; S.tTQ = S.t;
    S.Ip *= (1 + C.IPSPIKE);                       // pic d'aplatissement du profil
    S.prFlash = C.KW * p.n19 * S.Te0 / C.TQTAU;    // flash radiatif ~ ΔW/τ_TQ (W)
    S.neSpike = C.NEINFLUX * p.n19;
  }

  // --- Bilan d'énergie / dents de scie / ELM
  if (S.tq) {
    S.Te0 += (C.TQTE - S.Te0) * dt / C.TQTAU;
  } else {
    var fdeg = Math.max(C.FDEGMIN, 1 - C.FDEGK * (S.w + S.w32) / C.A);
    var wth = C.KW * p.n19 * S.Te0;                                  // J
    S.Te0 += (p.pheat * 1e6 - wth / (tauE * fdeg)) / (C.KW * p.n19) * dt;
    S.sawT += dt;
    if (S.sawT > S.sawP && S.Te0 > C.SAWTE) { S.sawT = 0; S.Te0 *= C.SAWDROP; }
    var elmA = p.elm;
    if (elmA > 0 && !S.locked && S.Te0 > C.ELMTE) {
      S.elmT -= dt;
      if (S.elmT <= 0) {
        S.elmB = elmA;
        S.Te0 *= (1 - C.ELMDROP * elmA);
        var jit = 1 + C.ELMJIT * Math.max(-1, Math.min(1, g()));
        S.elmT = C.ELMP * jit;
      }
    }
  }
  S.elmB *= Math.exp(-dt / C.ELMTAU);
  S.elmPh += 2 * Math.PI * C.ELMFHZ * dt;

  // --- Quench de courant : L/R avec résistivité post-quench
  if (S.tq && !S.cq && S.Te0 < C.CQTE) { S.cq = true; S.tCQ = S.t; S.TeCQ = S.Te0; }
  if (S.cq) {
    S.TeCQ += (C.TECQ - S.TeCQ) * dt / 0.002;      // relaxe vers ~10 eV
    var rp = etaSpitzer(S.TeCQ, C.ZEFFCQ, 1) * 2 * Math.PI * C.R0 /
             (Math.PI * C.A * C.A * C.KAPPA);
    S.Ip += (-S.Ip * rp / C.LP) * dt;
    if (S.Ip < C.ENDIP * C.IP0) { S.Ip = C.ENDIP * C.IP0; S.ended = true; S.tEnd = S.t; }
  }
  S.prFlash *= Math.exp(-dt / C.TQTAU);

  // --- Rayonnement et densité (canaux mesurés)
  var prEq = C.PRADBASE * p.pheat * 1e6 + C.PRADW * (S.w + S.w32) / C.A;
  S.Prad += (prEq - S.Prad) * dt / 0.01;
  S.Ne += ((S.cq ? -(1 - C.NECQ) : 0) - S.Ne) * dt / C.NETAU;
  S.neSpike *= Math.exp(-dt / C.NESPT);

  S.rotPh += S.Om * dt;
  stepNoise(S, p, g, dt);

  if (S.t < S.tLock || S.tLock < 0) { S.phaseId = (S.w >= 0.02 * C.A) ? 1 : 0; }
  if (S.locked) S.phaseId = 2;
  if (S.tq) S.phaseId = 3;
  if (S.cq) S.phaseId = 4;
  if (S.ended) S.phaseId = 5;
  S.t += dt;
  return measures(S, p);
}

// Pas de temps recommandé : fin pendant les quenches, grossier sinon.
function chooseDt(S) { return (S.tq && !S.ended) ? 1e-5 : 1e-4; }
/*MODEL-V6-END*/
  return {
    C: C,
    mulberry32: mulberry32,
    makeGauss: makeGauss,
    etaSpitzer: etaSpitzer,
    tauE98: tauE98,
    betaP: betaP,
    newState: newState,
    stepModel: stepModel,
    chooseDt: chooseDt
  };
});
