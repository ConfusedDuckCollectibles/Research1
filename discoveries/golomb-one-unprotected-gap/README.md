# When one gap controls a real Golomb ruler

A Golomb ruler places marks so that every pair gives a different distance. Here the marks can slide on the real line, but the distances must keep their original order and stay at least one unit apart. How short can such a ruler become?

If all but one adjacent mark gap lie within the ruler's initial run of consecutive occupied distances, an exact formula gives the answer. No search over possible real marks is needed: keep the two blocks of marks intact and shift one block by a single integer amount. The result applies to every cyclic perfect difference set cut at its second-largest circular gap. Even after compression, that cut cannot beat the largest-gap cut of the same cyclic set. An exact smaller-cut example shows that this comparison does not extend to all cuts.

The exact replay checked 3,398 distinct Singer second cuts for prime-power orders through 49, found 30 strict shortenings, and certified uniqueness within those tested chambers. Universal uniqueness and the global real-versus-integer Golomb optimum are not settled here. The known largest-gap theorem is prior work, not a new claim.

- [Proof and scope](THEOREM.md)
- Replay: `node check.js` (standard JavaScript); compare with `EXPECTED.json`.
- [Independent audit](independent/AUDIT.md) and replay: `node independent/check.js`.

This is a mathematical research note, not peer reviewed.
