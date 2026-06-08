# Part XXVII — The Goldbach Complex Tensor Integral on the 6N Skeleton

*Volume III of the Arithmetic Geodynamics programme on the 6N skeleton.*

Binary Goldbach, written as a self-correlation of the prime field: with `H` the
prime indicator, the **Goldbach complex tensor integral**

> I(2E) = Σ_x H(x)·H(2E−x) = r(2E),

the count of prime pairs summing to 2E. Its Hardy–Littlewood normalization is the
Goldbach singular series

> 𝔖_G(2E) = 2Π₂ · ∏_{q|E, q>2} (q−1)/(q−2),  Π₂ = ∏_{p>2}(1 − 1/(p−1)²) ≈ 0.6602 (twin-prime constant).

### Result — a positive floor on the heuristic main term

**Theorem.** The product defining Π₂ converges absolutely (O(1/q²) tail), and
`𝔖_G(2E) ≥ 2Π₂ ≈ 1.3203 > 0` for every even 2E, with equality iff E is a power of
two.

### Numerical confirmation to 2E = 10⁸ (the Goldbach comet)

- **Raw count grows** (no constant floor): min r(2E) rises from 22 at 2E=1024 to
  1.5×10⁵ at 2^26 — like 2E/(ln 2E)².
- **Normalized count is floored** above zero: the powers of two pin
  r(2E)(ln 2E)²/(2E) onto the singular-series floor (~0.74 at 10⁸, drifting slowly
  toward Π₂ via the finite-size 1+2/ln correction), reproducing the **Part XIII**
  singular-series collapse at 10⁸ scale.

### Scope (stated without hedging)

The floor bounds the **singular series** (the expected main term), **not** the
actual count. The Hardy–Littlewood asymptotic is itself conjectural, and proving
r(2E) > 0 for all even 2E **is** binary Goldbach — which is **open**. Positivity of
the main term is not positivity of the count. Ternary Goldbach (Helfgott) and Chen's
p+P₂ theorem are not reproved here. **This is the HL heuristic main term plus a
numerical check, not a step toward a proof.**

## Layout

```
.
├── paper/    Chen_6N_Paper27.{tex,pdf} + figure
├── figures/  fig_goldbach.{pdf,png}
├── data/     goldbach_comet.csv  (2E, r(2E), normalized, is_power_of_two)
├── code/
│   ├── exp_goldbach.py        # autocorrelation r(2E) to 2E=1e8 (NumPy 2.x, multiprocessing)
│   ├── fig_goldbach_make.py   # comet + normalized-floor figure
│   └── verify_goldbach.py     # Pi2 convergence, S_G floor 2*Pi2, raw-grows/normalized-floored
├── CITATION.cff · .zenodo.json · LICENSE (MIT)
```

## Reproducing

```bash
pip install numpy matplotlib
python code/verify_goldbach.py   # Pi2, the 2*Pi2 floor, raw-vs-normalized split (fast)
python code/exp_goldbach.py      # sieve + autocorrelation to 2E=1e8 (~100MB RAM; writes gb_*.npy)
python code/fig_goldbach_make.py # comet + floor figure
```

Expected: Π₂ ≈ 0.6602; 𝔖_G ≥ 2Π₂ with equality at powers of two; raw count grows,
normalized floored.

## Scope

A new large-scale **confirmation** for this programme, **not** a new theorem about
Goldbach. It recovers Part XIII and the classical Hardy–Littlewood singular series,
and does not bear on the open binary Goldbach conjecture. Continues Part XXVI
(doi:10.5281/zenodo.20593654).

## License

MIT — see `LICENSE`.
