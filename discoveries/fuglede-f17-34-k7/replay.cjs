'use strict';
const { spawnSync } = require('node:child_process');
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const read = name => JSON.parse(fs.readFileSync(path.join(__dirname, name), 'utf8'));
function run(name, args = []) {
  const p = spawnSync(process.execPath, [path.join(__dirname, name), ...args], {
    encoding: 'utf8', maxBuffer: 10 * 1024 * 1024
  });
  if (p.error) throw p.error;
  if (p.status !== 0) throw Error(name + ': ' + p.stderr);
  return JSON.parse(p.stdout);
}
const suite = process.argv.slice(2);
if (suite.length && (suite.length !== 1 || suite[0] !== '--independent')) {
  throw Error('Use node replay.cjs [--independent]');
}
if (suite.length) {
  for (const stage of ['arc', 'low', 'top-catalogue', 'octic']) {
    const prefix = 'independent-' + stage;
    assert.deepStrictEqual(run(prefix + '.cjs'), read(prefix + '-expected.json'));
    console.log(prefix + ' replay matched');
  }
} else {
  const expected = read('expected.json');
  for (const stage of ['g', '6', '7', '8', '9']) {
    assert.deepStrictEqual(run('verify-k7.cjs', [stage]), {stage, out: expected.branches[stage]});
    console.log(stage + ' replay matched');
  }
}
