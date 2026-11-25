# Exact replay

Run from this directory. The scripts use Node.js, C++17, or Python 3, with standard libraries only. Names inside each directory match the checked transcripts.

```sh
set -eu
(cd 01-reductions &&
  node check52.js | diff - check52.expected &&
  node independent.js | diff - independent-expected.json)
(cd 02-layers &&
  g++ -O3 -std=c++17 layer_pair.cpp -o /tmp/f13-layer-52 && /tmp/f13-layer-52 | diff - layer_pair.out &&
  g++ -O3 -std=c++17 independent.cpp -o /tmp/f13-layer-52-independent && /tmp/f13-layer-52-independent | diff - independent.out &&
  g++ -O3 -std=c++17 audit_meet_in_middle.cpp -o /tmp/f13-layer-52-audit && /tmp/f13-layer-52-audit | diff - audit_meet_in_middle.out &&
  python3 check_minors.py | diff - check_minors.out)
(cd 03-nonlayer &&
  node check52-nonlayer.js | diff - check52-nonlayer.expected &&
  { node independent-planar.js; node independent-palette.js; } | diff - independent-expected.txt)
(cd 04-both-k9 &&
  node check-both-k9.js | diff - check-both-k9.expected &&
  node independent.js | diff - independent-expected.json)
(cd 05-endgame &&
  node check-endgame.js | diff - check-endgame.expected &&
  node independent.js | diff - independent-expected.json)
```

The optional `node 04-both-k9/check-both-k9.js --cert` prints all 432 normalized left-annihilator witnesses. For each anchor, the checked certificate kills the corresponding coefficient columns and evaluates to one on the required right-hand side. Independent scripts do not need these supplied certificates.
