# A proof technique that carries from 13 to other prime grids

On the 13 × 13 × 13 wraparound grid, our earlier proof needed several arguments specific to 13. This result gives one step that works for **any prime side length greater than 3**.

**The general step.** If a non-tiling-sized shape passed Fuglede’s wave test on such a grid, its coordinates would have to satisfy a family of exact arithmetic cancellation rules through degree three. These rules apply to both the shape and its wave frequencies. The proof uses geometry rather than a grid search, and the underlying geometric bound is sharp.

**What changed for 17.** Classical obstructions leave seven candidate sizes: 34, 68, 102, 136, 153, 170 and 204. Two separate exact computations narrow the first one: a hypothetical 34-point example cannot put three of its points on a line or eight in a plane. This note does not eliminate any of the seven. A [subsequent proof excludes size 204](../fuglede-f17-204/), leaving six sizes in the checked ledger.

- [General, computer-free geometric proof](CUBIC-MOMENTS.md)
- [Exact 17-grid case breakdown and coverage](CHECKPOINT.md)
- Independent audits: [initial gates](INDEPENDENT-34-AUDIT.md) and [final tested gate](INDEPENDENT-K8-AUDIT.md)

The code needs Node.js only, with no input downloads. The transcripts are JSON; the first script prints five JSON objects, compared with the `first` array in `expected.json`.

```sh
node verify-cubic.cjs
node verify-f17-first.cjs
node verify-f17-k8.cjs
node independent-34.cjs
node independent-k8.cjs
```

The main theorem for prime grids remains open. The work is AI-assisted, independently checked within the project, and not externally peer reviewed.
