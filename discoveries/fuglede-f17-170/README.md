# No 170-point spectral set on the 17 × 17 × 17 grid

A 170-block shape can't tile this wraparound grid because 170 doesn't divide 4,913. It also can't pass Fuglede’s wave test. This is the second complete size ruled out on the 17-grid, after [204 points](../fuglede-f17-204/).

The proof first forces any candidate shape and its partner frequencies into ten full layers. In a family of slices through those layers, each ten-point profile must use at most nine slopes. Two exact searches agree: of the 7,125 symmetry-reduced combinations of positions and allowed slopes, only nine admit any profile, and none admits more than two distinct profiles after removing a common height. That restriction is too strong to hold in two different directions.

- [Proof and search coverage](PROOF.md)
- [Self-contained four-direction shadow lemma](FOUR-SPECIAL-119.md)
- [Independent proof and search audit](INDEPENDENT-AUDIT.md)

Run the standard-library Node.js replay from this directory:

```sh
node replay.cjs
```

The work is AI-assisted and checked within the project. It has not been externally peer reviewed. **Five candidate sizes remain: 34, 68, 102, 136 and 153.**
