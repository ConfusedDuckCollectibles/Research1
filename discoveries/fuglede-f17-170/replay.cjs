'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { spawnSync } = require('node:child_process');
function read(name) {
  return JSON.parse(fs.readFileSync(path.join(__dirname, name), 'utf8'));
}
function run(name, args = []) {
  const p = spawnSync(process.execPath, [path.join(__dirname, name), ...args], {
    encoding: 'utf8', maxBuffer: 10 * 1024 * 1024
  });
  if (p.error) throw p.error;
  if (p.status !== 0) throw Error(name + ': ' + p.stderr);
  return JSON.parse(p.stdout);
}
for (const [stage, expected] of [
  ['verify-170-atlas.cjs', 'expected-atlas.json'],
  ...['atlas','planar17','shadow'].map(name => ['independent-' + name + '.cjs', 'independent-' + name + '-expected.json'])
]) {
  assert.deepStrictEqual(run(stage), read(expected));
  console.log(stage + ' matched');
}
for (const [n, profiles, variance, constant, boolean] of [
  [7, 8, 96, 0, 0], [8, 11, 7296, 24, 24], [9, 11, 7296, 24, 24]
]) {
  const counts = [0, ...Array(16).fill(profiles)];
  assert.deepStrictEqual(run('verify-four-special.cjs', [String(n)]), {counts, variance, constant, boolean});
  console.log('four-special n=' + n + ' matched');
}
