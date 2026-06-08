#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Goldbach interference comparator: r(2E) = #{p+q=2E, p,q prime} = autocorrelation
of the prime field, to 2E ~ 1e8. NumPy 2.x clean; multiprocessing over sampled 2E."""
import numpy as np, math, os, time
from multiprocessing import Pool

XMAX = 100_000_000
t0 = time.time()
def tlog(m): print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

# global prime sieve (inherited by fork workers, copy-on-write)
tlog(f"building prime sieve to {XMAX:,}")
SIEVE = np.ones(XMAX+1, dtype=np.bool_); SIEVE[:2] = False
for i in range(2, int(math.isqrt(XMAX))+1):
    if SIEVE[i]:
        SIEVE[i*i::i] = False
tlog("sieve done")

def goldbach_count(twoE):
    """unordered # of prime pairs p<=q with p+q = twoE"""
    E = twoE // 2
    a = SIEVE[2:E+1]                 # p = 2..E
    b = SIEVE[twoE-2:E-1:-1]         # twoE-p for p = 2..E  (length E-1)
    n = min(len(a), len(b))
    return int(np.count_nonzero(a[:n] & b[:n]))

def work(twoE):
    g = goldbach_count(twoE)
    R = g * math.log(twoE)**2 / twoE     # normalized (Hardy-Littlewood singular-series scale)
    return (twoE, g, R)

if __name__ == "__main__":
    # sample set: log-spaced comet + powers of 2 (singular-series floor) + 2*prime
    rng = np.random.default_rng(0)
    grid = np.unique((np.logspace(4, 8, 1400).astype(np.int64)//2)*2)
    grid = grid[(grid>=1000) & (grid<=XMAX)]
    pow2 = np.array([2**k for k in range(10,27) if 2**k<=XMAX], dtype=np.int64)   # 2E=2^k -> E=2^{k-1}, empty odd product
    # 2E = 2*p (E prime): pick primes near a log spread
    pr = np.array([p for p in (10007,100003,1000003,10000019,30000001,50000017) if SIEVE[p]], dtype=np.int64)
    twoP = 2*pr
    samples = np.unique(np.concatenate([grid, pow2, twoP]))
    samples = samples[samples<=XMAX]
    tlog(f"{len(samples)} sampled even centres; computing on {os.cpu_count()} cores")
    with Pool() as pool:
        res = pool.map(work, samples.tolist(), chunksize=8)
    res.sort()
    arr = np.array(res, dtype=np.float64)   # cols: 2E, g, R
    tlog("done computing")

    # floor analysis
    twoPi2 = 2*0.6601618158
    is_pow2 = np.array([ (int(x)&(int(x)-1))==0 for x in arr[:,0] ])
    Rmin_overall = arr[:,2].min()
    print(f"\n2*Pi2 (twin-prime const x2) = {twoPi2:.4f}")
    print(f"min normalized R over all samples         = {Rmin_overall:.4f}")
    print(f"min normalized R on powers of two (floor) = {arr[is_pow2,2].min():.4f}  "
          f"(mean {arr[is_pow2,2].mean():.4f})")
    print(f"raw count g: min={arr[:,1].min():.0f} at 2E={arr[arr[:,1].argmin(),0]:.0f}; "
          f"max={arr[:,1].max():.0f} at 2E={arr[arr[:,1].argmax(),0]:.0f}")
    # raw minimum grows? show g for powers of 2 across scale
    print("\n raw count on powers of two (shows the RAW infimum GROWS):")
    for x,g,R in arr[is_pow2][:, :]:
        print(f"   2E=2^{int(round(math.log2(x))):<2d}={int(x):>12,}  g={int(g):>7,}  R={R:.4f}")
    np.save("gb_data.npy", arr); np.save("gb_pow2.npy", is_pow2)
    tlog("saved gb_data.npy")
