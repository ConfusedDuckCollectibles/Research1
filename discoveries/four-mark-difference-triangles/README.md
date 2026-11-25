# Four-mark difference triangle sets with minimum scope and minimum total span

We give exact constructions with **16, 17, 18, 20 and 21** four-mark rulers. In each construction, every difference between two marks of the same ruler occurs only once across all rulers, and the differences are precisely the integers from 1 to 6n. The sum of the ruler lengths also attains the lower bound 5n²+n.

| Number of rulers | Minimum scope | Minimum total span |
| ---: | ---: | ---: |
| 16 | 96 | 1296 |
| 17 | 102 | 1462 |
| 18 | 108 | 1638 |
| 20 | 120 | 2020 |
| 21 | 126 | 2226 |

Applying the **existing Wild product** to the 21-ruler construction gives an infinite sequence of perfect sets whose quadratic total-span coefficient improves the construction in *Higher-Order Staircase Codes*, Proposition 8, from 5+1/91 to **5+1/127**. This is a comparison with that specific construction, not a universal optimality or priority claim. The general simultaneous problem remains open here; no claim is made for the missing size 19.

Run the following to check the five certificates and a product with 5,334 rulers using only the standard library:

    python3 replay.py

See [the proof](PROOF.md) for the exact-cover reduction, product and primary literature. The file search.cpp is the original bounded discovery program, not part of the trusted base for the witnesses. This note and code were AI-assisted and have not been peer reviewed.
