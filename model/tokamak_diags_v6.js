/* Couche diagnostics v6 (palier B) — transforme l'état physique 0D en signaux
   de salle de contrôle réalistes, comme les bases d'entraînement ML
   professionnelles (APODIS, DPRF, FRNN) :

   - Réseau poloïdal de bobines de Mirnov (dB/dt, T/s) : chaque bobine voit
     cos(m·θ_j − φ_mode) — la structure m est DANS les données, la
     décomposition en nombres de mode est à la charge de l'analyse, comme sur
     une vraie machine. Signal ∝ Ω : muet après verrouillage.
   - Boucles à selle toroïdales (δB_r, T) : champ STATIQUE ∝ w² — c'est le
     « locked mode detector » réel : il continue de voir l'îlot verrouillé
     que les Mirnov ne voient plus.
   - Canaux ECE radiaux (keV) : profil Te(r), dents de scie INVERSÉES autour
     du rayon d'inversion (le cœur chute, l'extérieur reçoit le pulse de
     chaleur), aplatissement local au passage de l'îlot sur q=2.
   - Bolométrie (MW) : puissance rayonnée filtrée par le temps de réponse du
     détecteur (le flash du quench est lissé et retardé).
   - Interférométrie (1e19 m⁻²) : densité intégrée sur corde centrale.
   - Courant plasma (MA).

   Chaque canal porte son propre bruit AR(1) (correct en dt d'échantillonnage).
   Tout est déterministe par seed : les tirages passent par le gauss fourni.
   Double export : window.TokamakDiagsV6 / module.exports. */
(function (root, factory) {
  if (typeof module === 'object' && typeof module.exports === 'object') {
    module.exports = factory();
  } else {
    root.TokamakDiagsV6 = factory();
  }
})(typeof self !== 'undefined' ? self : this, function () {
/*DIAGS-V6-BEGIN*/
var DCFG = {
  NCOILS: 8,            // bobines de Mirnov, réparties poloïdalement
  BWALL: 1.15,          // rayon des capteurs (fraction de a)
  NSADDLE: 4,           // boucles à selle, réparties toroïdalement (n=1)
  CSAD: 0.25,           // calibration δB_r paroi (→ ~2 mT à w = 10 cm, JET-like)
  ECER: [0.10, 0.25, 0.40, 0.55, 0.70, 0.85], // rayons ECE (fractions de a)
  ECEALPHA: 1.2,        // forme du profil Te(r) = Te0·(1−(r/a)²)^α
  RINV: 0.40,           // rayon d'inversion des dents de scie (fraction de a)
  ECETAU: 0.03,         // relaxation des canaux ECE vers le profil (s)
  HEATPULSE: 0.25,      // fraction du ΔTe0 de crash convertie en pulse externe
  HEATTAU: 0.02,        // décroissance du pulse de chaleur (s)
  FLATMAX: 0.5,         // aplatissement max du canal ECE dans l'îlot
  TAUBOLO: 0.004,       // temps de réponse du bolomètre (s)
  CHORD: 0.8,           // facteur de corde de l'interféromètre
  ELMBB: 0.8,           // amplitude large-bande des ELM sur les bobines (T/s)
  // Écarts-types de bruit à noise = 1.5 % (mis à l'échelle par P.noise/1.5)
  SIGCOIL: 0.05,        // T/s
  SIGSAD: 2e-5,         // T
  SIGECE: 0.015,        // relatif au Te local (+ plancher 5 eV)
  SIGBOLO: 0.2,         // MW
  SIGNEL: 0.05,         // 1e19 m^-2
  SIGIP: 0.015,         // relatif à Ip0
  NOISETAU: 1.5e-3,     // corrélation AR(1) (s)
};

function profil(rA) { return Math.pow(Math.max(0, 1 - rA * rA), DCFG.ECEALPHA); }

// Crée l'état de la couche diagnostics pour un tir (S : état modèle v6 initial).
function newDiag(P, S) {
  var M = S.M, k;
  var D = {
    M: M,
    thetas: [], phis: [],
    ece: [], heat: [],
    bolo: 0, prevTe0: S.Te0,
    // états de bruit AR(1) par canal
    nCoil: [], nSad: [], nEce: [], nBolo: 0, nNel: 0, nIp: 0,
    geoDecay21: Math.pow(M.RS21 / DCFG.BWALL, 3),   // (r_s/b)^(m+1), m=2
    geoDecay32: Math.pow(M.RS32 / DCFG.BWALL, 4),   // m=3
  };
  for (k = 0; k < DCFG.NCOILS; k++) { D.thetas.push(2 * Math.PI * k / DCFG.NCOILS); D.nCoil.push(0); }
  for (k = 0; k < DCFG.NSADDLE; k++) { D.phis.push(2 * Math.PI * k / DCFG.NSADDLE); D.nSad.push(0); }
  for (k = 0; k < DCFG.ECER.length; k++) {
    D.ece.push(S.Te0 * profil(DCFG.ECER[k]));
    D.heat.push(0);
    D.nEce.push(0);
  }
  return D;
}

function channelNames() {
  var names = ['t_s'], k;
  for (k = 0; k < DCFG.NCOILS; k++) names.push('mir_p' + String(Math.round(360 * k / DCFG.NCOILS)).padStart(3, '0'));
  for (k = 0; k < DCFG.NSADDLE; k++) names.push('sad_t' + String(Math.round(360 * k / DCFG.NSADDLE)).padStart(3, '0'));
  for (k = 0; k < DCFG.ECER.length; k++) names.push('ece_r' + String(Math.round(100 * DCFG.ECER[k])).padStart(3, '0'));
  names.push('bolo_MW', 'nel_1e19m2', 'ip_MA');
  return names;
}

// Échantillonne tous les canaux. dts = période d'échantillonnage (s).
// p = paramètres résolus du tir (defP du modèle). Retourne un tableau de
// valeurs alignées sur channelNames().
function sampleDiag(D, S, p, g, dts) {
  var M = D.M, k, out = [S.t];
  var rho = Math.exp(-dts / DCFG.NOISETAU);
  var qf = Math.sqrt(1 - rho * rho) * (p.noise / 1.5);
  var wa = S.w / M.A, w32a = S.w32 / M.A;
  var alive = !S.ended;

  // --- Amplitudes de champ à la paroi
  var b21 = alive ? DCFG.CSAD * M.B0 * wa * wa * D.geoDecay21 : 0;        // T
  var b32 = alive ? DCFG.CSAD * M.B0 * w32a * w32a * D.geoDecay32 : 0;
  var ph21 = S.rotPh, ph32 = 1.5 * S.rotPh + 1.1;

  // --- Bobines de Mirnov : dB/dt ∝ Ω, structure poloïdale m
  for (k = 0; k < DCFG.NCOILS; k++) {
    var th = D.thetas[k];
    var sig = b21 * S.Om * Math.sin(2 * th - ph21)
            + b32 * S.Om * Math.sin(3 * th - ph32);
    if (S.elmB > 1e-3) sig += DCFG.ELMBB * S.elmB * g();  // bouffée large-bande
    D.nCoil[k] = rho * D.nCoil[k] + qf * DCFG.SIGCOIL * g();
    out.push(sig + D.nCoil[k]);
  }

  // --- Boucles à selle : δB_r statique, structure toroïdale n=1 (le 2/1
  //     domine ; le 3/2 est n=2, invisible sur une différence n=1 parfaite —
  //     on ajoute sa fuite géométrique à 10 %).
  for (k = 0; k < DCFG.NSADDLE; k++) {
    var phi = D.phis[k];
    var br = b21 * Math.cos(phi - ph21) + 0.1 * b32 * Math.cos(2 * phi - ph32);
    D.nSad[k] = rho * D.nSad[k] + qf * DCFG.SIGSAD * g();
    out.push(br + D.nSad[k]);
  }

  // --- ECE : crash de dents de scie inversé autour de RINV
  var crash = !S.tq && S.Te0 < D.prevTe0 * 0.93;
  var dTe0 = D.prevTe0 - S.Te0;
  for (k = 0; k < DCFG.ECER.length; k++) {
    var rA = DCFG.ECER[k];
    var cible = S.Te0 * profil(rA);
    // aplatissement local au passage de l'îlot sur q=2
    var demiW = S.w / (2 * M.A);
    if (demiW > 0.01) {
      var recouvre = Math.max(0, 1 - Math.abs(rA - M.RS21) / demiW);
      if (recouvre > 0) {
        var teIlot = S.Te0 * profil(M.RS21);
        cible = cible + (teIlot - cible) * DCFG.FLATMAX * recouvre;
      }
    }
    if (crash) {
      if (rA > DCFG.RINV) D.heat[k] += DCFG.HEATPULSE * dTe0 * (rA - DCFG.RINV) / (1 - DCFG.RINV);
      else D.ece[k] = cible; // le cœur suit le crash instantanément
    }
    var tau = S.tq ? 3e-4 : DCFG.ECETAU;
    D.ece[k] += (cible - D.ece[k]) * Math.min(1, dts / tau);
    D.heat[k] *= Math.exp(-dts / DCFG.HEATTAU);
    D.nEce[k] = rho * D.nEce[k] + qf * (DCFG.SIGECE * Math.max(cible, 0.3) + 0.005) * g();
    out.push(D.ece[k] + D.heat[k] + D.nEce[k]);
  }
  D.prevTe0 = S.Te0;

  // --- Bolométrie : filtre du premier ordre sur la puissance rayonnée totale
  var pradW = S.Prad + S.prFlash;
  D.bolo += (pradW - D.bolo) * Math.min(1, dts / DCFG.TAUBOLO);
  D.nBolo = rho * D.nBolo + qf * DCFG.SIGBOLO * g();
  out.push(D.bolo / 1e6 + D.nBolo);

  // --- Interférométrie : densité intégrée sur corde centrale
  var ne19 = p.n19 * (S.Ne + 1) + S.neSpike;
  D.nNel = rho * D.nNel + qf * DCFG.SIGNEL * g();
  out.push(ne19 * 2 * M.A * DCFG.CHORD + D.nNel);

  // --- Courant plasma
  D.nIp = rho * D.nIp + qf * DCFG.SIGIP * (M.IP0 / 1e6) * g();
  out.push(S.Ip / 1e6 + D.nIp);

  return out;
}
/*DIAGS-V6-END*/
  return {
    DCFG: DCFG,
    profil: profil,
    newDiag: newDiag,
    channelNames: channelNames,
    sampleDiag: sampleDiag
  };
});
