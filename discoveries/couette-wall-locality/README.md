# What the walls change in a mixing flow

Picture a channel where fluid layers slide past one another and stretch a streak of dye until diffusion smooths it out. How much of the most persistent streak can survive? In this linear model, a wide channel can be analyzed one wall at a time, even if the flow speeds up, slows down, or reverses.

The proof gives a sharp worst-case error for replacing the full channel with the more persistent of its two single-wall models. It also identifies the exact threshold for a genuinely wall-localized pattern. Above that threshold the approximation gets exponentially better with channel width. In one exactly checked reflecting-wall example, reversing the flow leaves **over four times as much worst-case residual** as the same model without walls predicts. This compares operator norms, not a particular streak, an experimental signal, or turbulent flow.

- [Proof](PROOF.md) and [closest prior studies](PRIOR.md)
- [Exact replay](REPLAY.md): run `python3 verify.py` (Python standard library only).
- [Independent audit](audit/AUDIT.md), including a sharp correction to the reflecting-wall kernel bound. Its separate numerical-free certificate is `python3 audit/sharp_constant.py`.
