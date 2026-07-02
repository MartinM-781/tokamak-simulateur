// Tests du générateur de datasets (scripts/generate.js) : cohérence du
// manifest, conformité du format CSV à l'export web, et déterminisme complet.
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const { execFileSync } = require('node:child_process');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');

const GEN = path.join(__dirname, '..', 'scripts', 'generate.js');

function genere(dir, args) {
  return execFileSync(process.execPath, [GEN, '--out', dir, ...args], { encoding: 'utf8' });
}

test('générateur : manifest cohérent, fichiers présents, format CSV conforme', () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_gen_'));
  try {
    genere(dir, ['--shots', '6', '--disrupt-ratio', '0.5', '--seed-base', '777']);
    const man = JSON.parse(fs.readFileSync(path.join(dir, 'manifest.json'), 'utf8'));
    assert.equal(man.tirs.length, 6);
    assert.equal(man.tirs.filter(t => t.label === 'disruptif').length, 3);
    assert.equal(man.dt_ms, 0.05);

    for (const tir of man.tirs) {
      const csv = fs.readFileSync(path.join(dir, tir.fichier), 'utf8');
      const lignes = csv.trim().split('\n');
      // En-têtes identiques à l'export de web/tokamak_v5.html
      assert.equal(lignes[0], '# tokamak_v5 synthetic shot');
      assert.match(lignes[1], /^# seed=\d+ delta0=-?[\d.]+ f0_khz=[\d.]+ cb=[\d.]+ noise_pct=[\d.]+$/);
      assert.match(lignes[2], /^# t_lock_ms=(NA|[\d.]+) t_tq_ms=(NA|[\d.]+) t_cq_ms=(NA|[\d.]+)$/);
      assert.equal(lignes[3], 't_ms,mirnov,te,ip,prad,ne,w21,w32,k_chirikov,f_khz,locked,phase');
      assert.match(lignes[4], /^\d+\.\d{2},-?\d+\.\d{5},-?\d+\.\d{5},-?\d+\.\d{5},-?\d+\.\d{5},-?\d+\.\d{5},\d+\.\d{5},\d+\.\d{5},\d+\.\d{4},-?\d+\.\d{4},[01],[0-5]$/);
      assert.equal(lignes.length - 4, tir.n_lignes_csv);
      // Le label visé correspond à l'issue réelle, elle-même cohérente avec t_tq
      assert.equal(tir.disrupte, tir.t_tq_ms !== null);
      assert.equal(tir.disrupte, tir.label === 'disruptif',
        `tir ${tir.seed} : label=${tir.label} mais disrupte=${tir.disrupte}`);
      // Un tir disruptif s'arrête avant TMAX, un tir sain couvre les 2500 ms
      if (tir.label === 'disruptif') {
        assert.ok(tir.t_end_ms !== null && tir.t_end_ms < 2500);
        // Ordre structurel seulement : le CQ suppose le TQ, la fin suppose le
        // CQ. L'ordre tLock < tTQ, lui, n'est PAS universel : à d0 et cb
        // faibles l'îlot peut atteindre K=1 avant le verrouillage (observé
        // sur ~2 % des tirs disruptifs, inversion ≈ 1,6 ms).
        assert.ok(tir.t_tq_ms < tir.t_cq_ms && tir.t_cq_ms < tir.t_end_ms);
      } else {
        assert.equal(tir.t_end_ms, null);
        assert.ok(tir.n_pas >= 2500 / 0.05 - 1);
      }
    }
  } finally {
    fs.rmSync(dir, { recursive: true, force: true });
  }
});

test('générateur : déterminisme — mêmes arguments ⇒ sorties identiques octet pour octet', () => {
  const d1 = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_gen_a_'));
  const d2 = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_gen_b_'));
  try {
    const args = ['--shots', '4', '--disrupt-ratio', '0.5', '--seed-base', '55'];
    genere(d1, args);
    genere(d2, args);
    for (const f of fs.readdirSync(d1).sort()) {
      const a = fs.readFileSync(path.join(d1, f), 'utf8');
      const b = fs.readFileSync(path.join(d2, f), 'utf8');
      assert.equal(a, b, `divergence sur ${f}`);
    }
    assert.equal(fs.readdirSync(d1).length, fs.readdirSync(d2).length);
  } finally {
    fs.rmSync(d1, { recursive: true, force: true });
    fs.rmSync(d2, { recursive: true, force: true });
  }
});

test('générateur : seed-base différent ⇒ dataset différent', () => {
  const d1 = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_gen_c_'));
  const d2 = fs.mkdtempSync(path.join(os.tmpdir(), 'tok_gen_d_'));
  try {
    genere(d1, ['--shots', '2', '--disrupt-ratio', '0.5', '--seed-base', '55']);
    genere(d2, ['--shots', '2', '--disrupt-ratio', '0.5', '--seed-base', '56']);
    const a = fs.readFileSync(path.join(d1, 'manifest.json'), 'utf8');
    const b = fs.readFileSync(path.join(d2, 'manifest.json'), 'utf8');
    assert.notEqual(a, b);
  } finally {
    fs.rmSync(d1, { recursive: true, force: true });
    fs.rmSync(d2, { recursive: true, force: true });
  }
});
