# Modèle v6 — 0D en unités SI, calibré sur la littérature

## Contrat épistémique

La v6 ne « valide » pas une théorie sur les vrais tokamaks : c'est un modèle
0D dont la chaîne causale (tearing → verrouillage → recouvrement → quench
thermique → quench de courant) est structurelle. Ce qu'elle apporte par
rapport à la v5 pédagogique :

1. **Unités réelles** (m, s, T, MA, keV, MW) sur une machine JET-like.
2. **Chaque terme est une expression publiée**, pas une constante d'ambiance :
   résistivité de Spitzer, temps résistif, scaling de confinement IPB98(y,2),
   décroissance L/R du courant avec inductance et résistivité post-quench.
3. **Les tests comparent le modèle à des chiffres expérimentaux publiés**
   (voir tests/model_v6.test.js) — le modèle devient falsifiable contre la
   littérature, pas contre lui-même.

Toute conclusion tirée de la v6 vaut *dans les hypothèses du modèle*
(0D, profils figés, tearing classique) et doit être énoncée comme telle.

## Machine de référence (préréglage JET-like, surchargeable)

| Grandeur | Valeur | Commentaire |
|---|---|---|
| R₀ / a | 3.0 m / 1.0 m | rapport d'aspect 3, section circulaire (κ=1) |
| B₀ | 3.0 T | champ toroïdal |
| Ip₀ | 2.5 MA | courant plasma de plat |
| n̄e | 3×10¹⁹ m⁻³ | paramètre de tir |
| P_chauffage | 10 MW | paramètre de tir |
| Z_eff / lnΛ | 1.7 / 17 | résistivité |
| r_s (q=2) / r_s' (q=1.5) | 0.75 a / 0.60 a | surfaces rationnelles (profil q figé) |
| τ_w (paroi) | 5 ms | temps résistif de la chambre |

## Équations

**Résistivité (Spitzer ∥, correction néoclassique)** :
η = NEO · 1.65×10⁻⁹ · lnΛ / Te[keV]^{3/2} Ω·m, NEO = 2.5 (piégeage à mi-rayon).
Te(r_s) = 0.4·Te0 (facteur de profil).

**Îlot 2/1 (équation de Rutherford, forme classique)** :
dw/dt = 1.22·(η/μ₀)·Δ′_eff, avec Δ′_eff·a = d0a·(1 − w/w_sat)
+ boost de verrouillage (perte d'écrantage rotationnel, +1.5 après lock).
Terme bootstrap NTM disponible (P.bs, défaut 0 — module de recherche) :
Δ′_bs·a = bs·√ε_s·β_p·(w/a)/((w/a)² + (w_d/a)²), w_d = 0.02 a.
Couplage 2/1 → 3/2 hérité de la v4/v5, en Δ′ dimensionné.

**Rotation et verrouillage (couple de paroi résistive, forme de Fitzpatrick)** :
dΩ/dt = (Ω₀−Ω)/τ_mom − (c_w·cwf/τ_w)·(w/a)⁴·Ω/(1 + Ω·τ_w),
τ_mom = τ_E. Verrouillage quand Ω < 0.05·Ω₀ (bifurcation par emballement :
le freinage ∝ w⁴ croît avec l'îlot et diverge quand Ω chute).
c_w est calibré pour verrouiller à w ≈ 2–5 cm (ordre des seuils observés).

**Bilan d'énergie** : W_th = k_W·n̄₁₉·Te0, k_W = 3·e·10¹⁹·V/2.5 (profil piqué).
dTe0/dt = (P_chauff − W_th/(τ_E·f_deg))/(k_W·n̄₁₉), τ_E = IPB98(y,2) :
τ_E = 0.0562·Ip^0.93·B^0.15·n̄^0.41·P^−0.69·R^1.97·ε^0.58·κ^0.78·M^0.19
(Ip MA, P MW, n̄ 10¹⁹ m⁻³). f_deg = max(0.15, 1 − 1.5·(w+w₃₂)/a) :
dégradation du confinement par les îlots.

**Dents de scie** : période 120 ms (±10 % par tir), crash Te0 ×0.85 tant que
le cœur est chaud (Te0 > 2 keV ⇔ proxy q₀ < 1). **ELM** : période ~30 ms
(~33 Hz, type-I) jitterée, chute de W_th de ~2 %·elm, bouffée Mirnov + P_rad.

**Quench thermique** : déclenché par recouvrement K = (w+w₃₂)/(r_s−r_s') ≥ 1
(critère de Chirikov, v4). Effondrement de Te0 en τ_TQ = 0.3 ms vers 50 eV,
flash radiatif P_rad ≈ ΔW_th/τ_TQ (échelle GW).

**Quench de courant** : pic d'Ip de +5 % au TQ (aplatissement du profil),
puis dIp/dt = −Ip·R_p/L_p avec L_p = μ₀R₀(ln(8R₀/a) − 2 + ℓᵢ/2), ℓᵢ = 1,
R_p = η(Te_CQ, Z_eff=3)·2πR₀/(πa²κ), Te_CQ relaxant vers 10 eV.
Ordres de grandeur obtenus : τ_L/R ≈ 13 ms, quench 80→20 % normalisé à la
section ≈ 5–6 ms/m² — comparé au plancher ITER de 1.7 ms/m² et aux
distributions JET (tests).

**Mesures** (échantillonnées à 1 kHz, bruit AR(1) correct en dt :
ρ = exp(−dt/τ_n), τ_n = 1.5 ms) :
- `mirnov_Tps` : dB/dt paroi ≈ c_mir·B₀·(w/a)²·Ω·sin(φ) + composante 3/2 —
  ~qq T/s pour w ≈ 3 cm à qq kHz (ordre réaliste des bobines de Mirnov) ;
- `te0_keV` (ECE central), `ip_MA`, `prad_MW` (bolométrie), `ne_1e19`
  (interférométrie, bouffée d'influx au TQ).

## Références utilisées pour les cibles de tests

- ITER Physics Basis (Nucl. Fusion 39, 1999) : scaling IPB98(y,2) ;
  plancher du quench de courant normalisé 1.7 ms/m².
- de Vries et al., Nucl. Fusion 51 (2011) : statistiques des causes et
  précurseurs de disruptions à JET (durées de précurseur ~0.1–qq s).
- Sweeney et al., Nucl. Fusion 57 (2017) : base de modes verrouillés DIII-D
  (ordres de grandeur des seuils et délais lock → disruption).
- Wesson, « Tokamaks » : Spitzer, inductance, Rutherford, dents de scie.

Les plages encodées dans les tests sont volontairement larges (ordres de
grandeur et bandes publiées, pas des points) : elles échouent si la physique
dimensionnée dérive, pas si un détail numérique bouge.
