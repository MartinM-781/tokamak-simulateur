# Tokamak — interactive plasma physics & disruption prediction

[![CI](https://github.com/MartinM-781/tokamak-simulateur/actions/workflows/ci.yml/badge.svg)](https://github.com/MartinM-781/tokamak-simulateur/actions/workflows/ci.yml)
[![License MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Zero dependencies](https://img.shields.io/badge/dependencies-zero-brightgreen.svg)](package.json)
[![Live demo](https://img.shields.io/badge/▶-live%20demo-orange.svg)](https://martinm-781.github.io/tokamak-simulateur/web/index.html)

**English** | [Français](README.fr.md)

> Five interactive tokamak-physics simulators (self-contained HTML, zero dependencies),
> a deterministic generator of labelled synthetic disruption datasets, and a Python ML
> baseline for early disruption detection — all held together by physics regression tests.

**[▶ Open the live demo](https://martinm-781.github.io/tokamak-simulateur/web/index.html)** —
the five simulators, running directly in your browser. No install, no build, no server.

![Act V — a full disruptive shot: Mirnov precursor, mode locking (red zone), thermal quench then current quench](docs/captures/tokamak_v5.png)

A tokamak is the magnetic bottle that confines a fusion plasma at 100+ million degrees.
Its worst failure mode is a **disruption** — the plasma loses control in milliseconds and
dumps its energy into the wall. This repo is a hands-on way to *see* the physics that leads
there, and then to *measure* how early a detector can call it. Every simulator runs from a
single HTML file (double-click, `file://`, vanilla JS); the dataset generator and the ML
baseline reproduce, bit-for-bit from a seed, the same physics you watched in the browser.

## The five acts

Each page opens **straight in the browser** (double-click — no server, no install) and
follows the causal chain of a disruption:

| | |
|---|---|
| [![v1](docs/captures/tokamak_v1.png)](web/tokamak_v1.html) **Act I — Field lines & safety factor**: helical winding, rational surfaces, why q = 2 is a fault line. | [![v2](docs/captures/tokamak_v2.png)](web/tokamak_v2.html) **Act II — Particle orbits**: Boris pusher, magnetic mirror, banana orbits, vertical drift with no plasma current. |
| [![v3](docs/captures/tokamak_v3.png)](web/tokamak_v3.html) **Act III — MHD equilibrium**: ∇P = j×B, Shafranov shift, operational limits (q(a), Troyon, density). | [![v4](docs/captures/tokamak_v4.png)](web/tokamak_v4.html) **Act IV — Magnetic islands**: tearing mode, island chains, Chirikov overlap criterion K ≥ 1. |

**[Act V — Disruption](web/tokamak_v5.html)** (screenshot above): Rutherford equation, mode
locking by wall torque, thermal quench then current quench. Every shot produces labelled
noisy signals, exportable to CSV. Index page: [web/index.html](web/index.html).

## The ML bench — predicting the disruption before it happens

The v5 physics engine ([model/tokamak_model.js](model/tokamak_model.js)) runs as-is under
Node: entire datasets are generated from the CLI, with the exact same physics as the web
page, reproducible byte-for-byte from a seed (the AR(1) measurement noise is the only
stochasticity, and it is seeded).

```bash
# v5 — arbitrary units, 5 channels at 4 kHz (one CSV per shot + manifest)
node scripts/generate.js --shots 200 --out data/run01 --disrupt-ratio 0.5 --seed-base 1000

# v6 — SI units (JET-like machine), 0D at 1 kHz
node scripts/generate_v6.js --shots 200 --out data/run_v6 --disrupt-ratio 0.5 --seed-base 2000

# v6-diag — 22 raw control-room channels at 10 kHz (CSV.gz):
# 8 Mirnov coils, 4 saddle loops, 6 ECE radii, bolometry, interferometry…
node scripts/generate_v6_diags.js --shots 100 --out data/run_v6_diag02 --disrupt-ratio 0.5 --seed-base 3000
```

The Python baseline ([ml/](ml/)) extracts windowed features, calibrates robust z-score
detectors (median/MAD over healthy shots) and a logistic regression, then evaluates
**under realistic conditions**: split by shot, thresholds calibrated on the train set only,
and on a disruptive shot only the windows *before* the thermal quench count.

```bash
python ml/features.py data/run01 && python ml/baseline.py data/run01 && python ml/evaluate.py data/run01
# _v6 variants for the diagnostics datasets
```

### The headline result: the locked-mode blind spot

A disruption precursor starts as a *rotating* mode, highly visible on the Mirnov coils
(dB/dt ∝ W²·Ω). But once the island **locks** (Ω → 0), the probe goes silent **while the
danger keeps growing** — this is the historical blind spot of purely magnetic detection,
and the datasets reproduce it:

| Detector (v6-diag, pre-quench windows only) | Rotating mode | Locked mode |
|---|---|---|
| `z_mirnov` (Mirnov array) | 97 % | **0 %** |
| `z_saddle` (saddle loops, static δB_r) | 100 % | **100 %** |
| `z_multi` (multi-diagnostic) | 100 % | 100 % |

The answer to the blind spot is not alarm memory, it is **the right sensor**: the static
radial field of the locked island stays visible on the saddle loops (the equivalent of the
locked-mode detector on real machines) — and in v5, it is the radiative rise (P_rad) that
takes over. The multi-channel detector inherits this physical complementarity: 100 %
detection, 0 % false alarms on the test set, median warning ~2 s before the quench.

![Mirnov vs saddle loops during the locked mode](ml/resultats_v6/mode_verrouille.png)

Full results, protocol and curves: [ml/RESULTS.md](ml/RESULTS.md) (v5) and
[ml/RESULTS_V6.md](ml/RESULTS_V6.md) (v6-diag).

## The v6 physics — SI units, falsifiable against the literature

The v6 engine ([model/tokamak_model_v6.js](model/tokamak_model_v6.js)) is a 0D model of a
JET-like machine (R₀ = 3 m, B₀ = 3 T, Ip = 2.5 MA, overridable presets) in which **every
term is a published expression**: Spitzer resistivity, Rutherford equation (+ optional NTM
bootstrap term), Fitzpatrick resistive-wall torque, IPB98(y,2) confinement scaling, L/R
current decay. The tests compare its outputs against published experimental numbers (ITER
Physics Basis, de Vries 2011, Sweeney 2017, Wesson) — the model is falsifiable against the
literature, not against itself. Equations and epistemic contract:
[docs/V6_PHYSIQUE.md](docs/V6_PHYSIQUE.md).

## Tests & validation

```bash
npm test          # 55 physics regression tests (node:test, zero dependencies)
npm run valide    # v6 battery: dt-convergence, extreme-range anti-NaN sweep,
                  # metamorphic properties, physical sensor bounds
```

Guarantees held by the tests: byte-for-byte determinism from a seed, purely observational
noise (the dynamics are identical at σ = 0 and σ = 5 %), strict parity between the web page's
CSV export and the CLI generator's, physical monotonicities (tTQ decreases with instability,
Te max increases with heating…).

## Layout

```
web/       5 self-contained HTML pages (file://, vanilla JS, zero build)
model/     pure physics engines (v5, v6, v6 diagnostics) — dual browser/Node export
scripts/   CLI dataset generators + v6 validation battery
tests/     physics regression tests (node:test)
ml/        Python baseline (features, detectors, evaluation, results)
docs/      v6 physics, screenshots
data/      generated datasets (gitignored — regenerable identically from a seed)
```

## Quickstart

- **Simulators**: open `web/index.html` in a browser (or the
  [live demo](https://martinm-781.github.io/tokamak-simulateur/web/index.html)). That's it.
- **Datasets**: Node ≥ 18, zero npm dependencies.
- **ML**: Python ≥ 3.10 — `pip install numpy pandas scikit-learn matplotlib`.

## License

[MIT](LICENSE)
