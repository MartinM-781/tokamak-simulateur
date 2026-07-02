// Tests du générateur v6 : cohérence du manifest, format CSV SI, déterminisme.
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const { execFileSync } = require('node:child_process');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');

const GEN = path.join(__dirname, '..', 'scripts', 'generate_v6.js');

function genere(dir, args) {
  return execFileSync(process.execPath, [GEN, '--out', dir, ...args], { encoding: 'utf8' });
}

test('générateur v6 : manifest cohérent et format CSV SI conforme', () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_gen6_'));
  try {
    genere(dir, ['--shots', '4', '--disrupt-ratio', '0.5', '--seed-base', '888']);
    const man = JSON.parse(fs.readFileSync(path.join(dir, 'manifest.json'), 'utf8'));
    assert.equal(man.schema, 'tokamak_v6_run/1');
    assert.equal(man.tirs.length, 4);
    assert.equal(man.tirs.filter(t => t.label === 'disruptif').length, 2);
    assert.equal(man.machine.Ip0_MA, 2.5);

    for (const tir of man.tirs) {
      const csv = fs.readFileSync(path.join(dir, tir.fichier), 'utf8');
      const lignes = csv.trim().split('\n');
      assert.equal(lignes[0], '# tokamak_v6 synthetic shot (SI, JET-like)');
      assert.match(lignes[1], /^# seed=\d+ d0a=-?[\d.]+ f0_khz=[\d.]+ cwf=[\d.]+ pheat_mw=[\d.]+ n19=[\d.]+ noise_pct=[\d.]+$/);
      assert.match(lignes[2], /^# t_lock_s=(NA|[\d.]+) t_tq_s=(NA|[\d.]+) t_cq_s=(NA|[\d.]+) t_onset_s=(NA|[\d.]+)$/);
      assert.equal(lignes[3], 't_s,mirnov_Tps,te0_keV,ip_MA,prad_MW,ne_1e19,w_m,w32_m,k_chirikov,f_khz,locked,phase');
      assert.equal(lignes.length - 4, tir.n_lignes_csv);
      assert.equal(tir.disrupte, tir.t_tq_s !== null);
      assert.equal(tir.disrupte, tir.label === 'disruptif',
        `tir ${tir.seed} : label=${tir.label} mais disrupte=${tir.disrupte}`);
      if (tir.label === 'disruptif') {
        assert.ok(tir.t_tq_s < tir.t_cq_s && tir.t_cq_s < tir.t_end_s);
        assert.ok(tir.w_at_lock_m === null || (tir.w_at_lock_m > 0.005 && tir.w_at_lock_m < 0.15));
      } else {
        assert.equal(tir.t_end_s, null);
        // Tir sain : couvre les 10 s à 1 kHz.
        assert.ok(tir.n_lignes_csv >= 10e3 - 2, `n_lignes=${tir.n_lignes_csv}`);
      }
    }
  } finally {
    fs.rmSync(dir, { recursive: true, force: true });
  }
});

test('générateur v6 : déterminisme octet pour octet', () => {
  const d1 = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_gen6_a_'));
  const d2 = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_gen6_b_'));
  try {
    const args = ['--shots', '2', '--disrupt-ratio', '0.5', '--seed-base', '99'];
    genere(d1, args);
    genere(d2, args);
    for (const f of fs.readdirSync(d1).sort()) {
      assert.equal(
        fs.readFileSync(path.join(d1, f), 'utf8'),
        fs.readFileSync(path.join(d2, f), 'utf8'),
        `divergence sur ${f}`
      );
    }
  } finally {
    fs.rmSync(d1, { recursive: true, force: true });
    fs.rmSync(d2, { recursive: true, force: true });
  }
});
