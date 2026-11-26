# A tighter limit on 34-point shapes in the 17³ Fuglede grid

If a 34-point shape and its partner frequencies passed Fuglede’s wave test on a 17 × 17 × 17 wraparound grid, **no plane could hold seven of their points**. The same limit of six applies to the number of evenly spaced projection directions in any one pencil: at most six projective normals in a pencil can split either shape into seventeen parallel planes of two points each.

The plane bound follows from an earlier reduction and a known geometric theorem. For the pencil bound, we reduce the possible projection profiles to ten symmetry classes and check every possible maximal polynomial degree from six through nine using exact arithmetic. The degree-nine search kept ten candidates through the fifth-coefficient check; none passed the next coefficient. The [independent audit](INDEPENDENT-AUDIT.md) includes separate implementations for the planar check, profile inventories and the degree-six, seven and eight gates, and a fresh replay of the frozen degree-nine verifier. A published bound of [Ball](https://doi.org/10.1112/jlms/54.3.581) then gives at most 86 difference directions per member.

**What this does not prove:** Size 34 remains open. This tightens necessary conditions from the [earlier prime-grid checkpoint](../fuglede-prime-cubic-moments/); it does not eliminate a full candidate size. The [204-point](../fuglede-f17-204/) and [170-point](../fuglede-f17-170/) exclusions are separate results.

## Proof and exact replay

- [Proof, coverage, and per-degree checks](PROOF.md)
- [Independent audit and separate replay](INDEPENDENT-AUDIT.md)
- [Expected main results](expected.json)

With Node.js installed, run these from this directory. No packages or downloads are needed. Each full command takes a few minutes.

```sh
node replay.cjs
node replay.cjs --independent
```

Run one main stage with `node verify-k7.cjs g` (plane cap) or use `6`, `7`, `8`, or `9`. The work is AI-assisted, checked within the project, and not externally peer reviewed.
