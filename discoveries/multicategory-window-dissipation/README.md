# Hidden dissipation from short, multi-category records

A recorded process can look exactly the same played forward or backward, even when the underlying dynamics produce entropy. The question here is how much dissipation can be certified when the data retain only a category from each half of a short recording.

For one exact 15-state Markov model, the estimate that uses all first-exit categories certifies **at least 4.55 times** the earlier optimized lower bound from the **same** short-record summary. It certifies **at least 27.94%** of the model's true dissipation. A simple forward-versus-backward comparison and every shared binary recoding of the six outcomes give zero. In a related mathematical family, the ratio between the bounds grows without bound, although the absolute dissipation gets smaller in the chosen time units.

**When the full record is available, use it.** The earlier uncensored excursion method already certifies at least 88.14% in this model. The improvement applies to the short-record summary, not to all ways of measuring dissipation. This is a theoretical result from a constructed model, not an experiment or a ready-to-use estimator.

- [General classifier bound, exact model, and prior-art comparison](proof.md)
- [Stronger first-exit bound and independent audit](first_exit_tent_audit.md). It uses the same short-record summary, without recording exit times.
- [Exact arithmetic replay](verify.py): run `python3 verify.py` and `python3 verify_first_exit_tent.py` with the Python standard library.
