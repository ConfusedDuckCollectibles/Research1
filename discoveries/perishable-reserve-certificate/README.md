# What separate inventory totals can hide

Counting how much stock meets a compatibility requirement and how much lasts until a deadline doesn't tell you how many units do **both**. If they are the same units, separately reassuring totals can hide a shared constraint.

In a synthetic 20-request example, two inventories have exactly the same compatibility and expiry totals. One serves all **20** requests and the other can serve only **11**. The data also guarantee seven requests from either of two groups when considered alone, but cannot guarantee seven from both groups together.

The proof gives the exact worst case without listing every possible stock table. It also identifies a set of requests that can always be served once the true compatibility–expiry pairings are revealed and proves when multiple protected quotas can coexist. The mathematical report assumes nested eligibility, stock all present initially, fixed service times and offline matching. It is not guidance for medical allocation, blood handling or a real shortage. No compatible public holdout was found.

[Proof and scope](PROOF.md) · [Prior work and data limitations](SOURCES.md) · [Independent audit](independent/AUDIT.md)

    python3 -B replay.py | diff -u replay.txt -
    (cd independent && python3 -B check.py | diff -u REPLAY.txt -)

Standard-library Python 3.10+. No patient, clinical or third-party data are needed. AI-assisted research note, not peer reviewed; no priority claim.
