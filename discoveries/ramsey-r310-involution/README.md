# A symmetry ruled out for the 40-vertex Ramsey candidate

**Result.** A triangle-free graph on 40 vertices with no independent set of 10 cannot have a nonidentity involution that fixes more than 32 vertices. An involution is a symmetry that undoes itself when applied twice. The specific case with 34 fixed vertices was left open in the [October 2026 Unsolved Lab report](https://unsolvedlab.com/problems/ramsey-r310).

The proof looks at the 34 fixed vertices after the six moved vertices pair off. Three overlapping 27-vertex induced graphs must have at least 255 edges in total, but their shared degree budget permits at most 249. The argument is elementary once the [published minimum-edge data](https://arxiv.org/abs/1208.0501) are taken as input.

[Proof](PROOF.md) · [Independent audit](AUDIT.md) · `node independent.js | diff - expected.json`

This closes a particular symmetry class rather than the general problem. The known range remains [40 ≤ R(3,10) ≤ 41](https://arxiv.org/abs/2401.00392). The note was independently checked within this project; it has not been peer reviewed. No unrestricted priority claim is made.
