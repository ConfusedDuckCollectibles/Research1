"use strict";
const fs = require("node:fs");
const path = require("node:path");
const { spawnSync } = require("node:child_process");
const { isDeepStrictEqual } = require("node:util");

const all = [
  "initial-pencils", "eight-uniform", "seven-uniform", "planar-hexagon",
  "no-rich-family", "at-most-two-rich", "exactly-three-rich", "marked-quintic",
  "quadric-controls", "quartic-twisted-cubic", "genus-two",
  "catalecticant-controls"
];
const selected = process.argv.length > 2 ? process.argv.slice(2) : all;
if (selected.some(id => !all.includes(id))) {
  console.error("Choose gates from: " + all.join(", "));
  process.exit(2);
}
for (const id of selected) {
  const started = Date.now();
  const source = path.join(__dirname, "verify", id + ".cjs");
  const expected = JSON.parse(fs.readFileSync(path.join(__dirname, "expected", id + ".json"), "utf8"));
  const run = spawnSync(process.execPath, [source], { encoding: "utf8", maxBuffer: 64 * 1024 * 1024 });
  if (run.error || run.status !== 0) {
    process.stderr.write(run.stderr || "");
    throw run.error || new Error(id + ": exit status " + run.status + " (signal " + run.signal + ")");
  }
  const actual = run.stdout.trim().split(/\r?\n/).map(JSON.parse);
  if (!isDeepStrictEqual(actual, expected)) {
    const i = actual.findIndex((record, j) => !isDeepStrictEqual(record, expected[j]));
    const j = i < 0 ? Math.min(actual.length, expected.length) : i;
    throw new Error(id + ": records differ near " + j + " (actual count " + actual.length +
                    ", expected " + expected.length + ")");
  }
  console.log(id + " PASS (" + actual.length + " records, " + (Date.now() - started) + " ms)");
}
