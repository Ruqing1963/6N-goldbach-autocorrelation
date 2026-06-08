#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np, math, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
arr = np.load("gb_data.npy"); is_pow2 = np.load("gb_pow2.npy")
twoE, g, R = arr[:,0], arr[:,1], arr[:,2]
Pi2 = 0.6601618158
fig,(a,b)=plt.subplots(1,2,figsize=(13.5,5.2))
fig.suptitle("Goldbach autocorrelation r(2E)=#{p+q=2E} to 2E~1e8: "
             "raw count grows, normalized floor stays > 0",fontsize=12,fontweight="bold")
# A: comet (raw)
a.loglog(twoE[~is_pow2], g[~is_pow2], ".", ms=2.5, color="#9ab", alpha=.5, label="even centres 2E")
a.loglog(twoE[is_pow2], g[is_pow2], "o", ms=5, color="#b4341f", label="2E = 2^k (lower envelope)")
a.set_xlabel("even centre 2E"); a.set_ylabel("Goldbach count r(2E)")
a.set_title("(A) the comet: raw r(2E) GROWS (no constant floor)")
a.legend(fontsize=8.5); a.grid(alpha=.3,which="both")
# B: normalized R = r * ln^2(2E)/(2E)
b.semilogx(twoE[~is_pow2], R[~is_pow2], ".", ms=2.5, color="#9ab", alpha=.5, label="even centres 2E")
b.semilogx(twoE[is_pow2], R[is_pow2], "o-", ms=5, color="#b4341f", label="2E = 2^k (singular-series floor)")
b.axhline(Pi2, color="k", ls="--", lw=1, label=f"$\\Pi_2$ = {Pi2:.3f} (HL floor, slow limit)")
b.axhline(2*Pi2, color="0.5", ls=":", lw=1, label=f"$2\\Pi_2$ = {2*Pi2:.3f}")
b.set_xlabel("even centre 2E"); b.set_ylabel(r"normalized  r(2E)$\cdot\ln^2(2E)/(2E)$")
b.set_title("(B) normalized floor frozen > 0 (powers of 2), -> $\\Pi_2$ slowly")
b.set_ylim(0,2.6); b.legend(fontsize=8); b.grid(alpha=.3,which="both")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig("fig_goldbach.png",dpi=200); fig.savefig("fig_goldbach.pdf")
print("wrote fig_goldbach.{png,pdf}")
print(f"powers-of-2 normalized floor at largest scales: {R[is_pow2][-3:]}")
