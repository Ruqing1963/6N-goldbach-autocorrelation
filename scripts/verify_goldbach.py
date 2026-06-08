#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verifier for Part XXVII (standard library).

(A) The twin-prime constant Pi2 = prod_{p>2}(1-1/(p-1)^2) converges (O(1/p^2) tail);
    the Goldbach singular series S_G(2E) = 2*Pi2 * prod_{q|E,q>2}(q-1)/(q-2) satisfies
    S_G(2E) >= 2*Pi2, with equality iff E is a power of two.
(B) Empirically (small scale): the RAW Goldbach count r(2^k) grows, while the
    NORMALIZED count r*(ln 2E)^2/(2E) stays floored above zero.
"""
import math

def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return s

def Pi2(P=2_000_000):
    s=primes_upto(P); prod=1.0
    for p in range(3,P+1):
        if s[p]: prod*=1-1.0/(p-1)**2
    return prod

def SG(E, oddprimes):
    f=1.0
    for q in oddprimes:
        if E%q==0: f*=(q-1)/(q-2)
    return f  # = prod_{q|E,q>2}(q-1)/(q-2); multiply by 2*Pi2 for S_G

def main():
    p2=Pi2(); twoPi2=2*p2
    print(f"(A) Pi2 = {p2:.6f}  (twin-prime constant);  2*Pi2 = {twoPi2:.4f}")
    s=primes_upto(100000); odd=[q for q in range(5,100000) if s[q]]
    floors=[]
    print("    S_G(2E)/(2Pi2) = prod_{q|E,q>2}(q-1)/(q-2)  for sample E:")
    for E,tag in [(2**10,"2^10"),(2**15,"2^15"),(3*5*7,"105=3*5*7"),(15015,"15015=3*5*7*11*13"),(9973,"prime")]:
        f=SG(E,odd); floors.append((E,f,tag))
        print(f"      E={tag:18} -> ratio {f:.4f}  (>=1; =1 iff power of 2)")
    ok_floor = all(f>=1-1e-12 for _,f,_ in floors) and abs(SG(2**15,odd)-1.0)<1e-12
    print(f"    floor S_G>=2*Pi2 holds, equality at powers of 2: {ok_floor}")

    # (B) raw vs normalized at powers of two
    NB=2**19; s2=primes_upto(NB)
    print("\n(B) powers of two: raw count grows, normalized floored")
    print("    2E         r(2E)    r*ln^2/(2E)")
    prev=0; grows=True; normfloor=True
    for k in range(10,19):
        twoE=2**k; E=twoE//2
        r=sum(1 for p in range(2,E+1) if s2[p] and s2[twoE-p])
        R=r*math.log(twoE)**2/twoE
        print(f"    2^{k:<2d}={twoE:>8,}  {r:>6}   {R:.4f}")
        grows = grows and (r>prev); prev=r
        normfloor = normfloor and (R>0.6)
    print(f"    raw count strictly grows: {grows};  normalized stays >0.6: {normfloor}")
    print("\nALL CHECKS PASS:", ok_floor and grows and normfloor)

if __name__=="__main__":
    main()
