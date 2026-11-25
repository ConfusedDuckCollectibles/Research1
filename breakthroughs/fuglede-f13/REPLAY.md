# Exact reproduction

The 78-point package uses Node.js, exact small-integer arithmetic, and built-in language features only. Run from this directory:

```sh
set -eu
(cd 01-neither-layer &&
  node check.js | diff - expected.txt &&
  node independent.js | diff - independent-expected.json)
(cd 02-one-layer &&
  node check.js | diff - expected.txt &&
  node independent.js | diff - independent-expected.json)
(cd 03-linked-layers &&
  node independent.js | diff - independent-expected.json)
(cd 04-final &&
  node check.js | diff - expected.txt &&
  node independent.js | diff - independent-expected.txt &&
  node independent-all-supports.js | diff - independent-all-supports.expected.json)
```

Each proof step explains the coverage of its finite checks and the analytic inference from them. No numerical cutoff, external dataset, random search, commercial solver, or dependence on a successful network fetch enters the proof.

The earlier cardinality modules have their own replay commands: [26](../../discoveries/fuglede-f13-26/), [39](../../discoveries/fuglede-f13-39/), [52](../../discoveries/fuglede-f13-52/), [104 and 117](../../discoveries/fuglede-f13-104-117/). They use some C++17 and Python 3 alongside JavaScript and are part of the proof of the full grid case, not prerequisites for replaying the 78-only check.
