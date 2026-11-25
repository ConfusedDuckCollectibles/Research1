# Two equally accurate sleep summaries can miss different awakenings

Knowing how many half-minute chunks a sleep tracker labeled correctly can still leave you unable to tell whether it found the same stretches of wakefulness.

In an eight-hour synthetic example, two predictions agree with the same reference **92.5%** of the time. They have identical distributions of sleep and wake run lengths and identical joint counts of every short label sequence through four epochs. Yet under the specified overlap rule one identifies **0 of 8** reference wake bouts and the other identifies **4 of 8**.

A separate construction has both channels reporting the same wake minutes and bout counts, and has **99.17%** epoch agreement in both examples, with **2 of 12** versus **8 of 12** matched bouts. The larger gap depends on the prescribed half-overlap scoring rule. The two cases make different assumptions about what is held fixed.

For this event-level precision and recall, the useful minimum report is simple: how many eligible bouts each track contains and how many match in time. These totals alone don't report uncertainty across people. No compatible public aggregate holdout was found. The result does not evaluate a device or give a clinical classification.

[Proof and assumptions](PROOF.md) · [Primary literature](SOURCES.md) · [Independent audit](independent/AUDIT.md)

The original and separate exact replay scripts require only Python 3.10 or later. No data download is needed:

    python3 -B replay.py | diff -u replay.txt -
    (cd independent && python3 -B check.py | diff -u REPLAY.txt -)

AI-assisted research note, not peer reviewed; no priority claim.
