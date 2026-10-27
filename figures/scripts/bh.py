import numpy as np
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); r=p[o]*n/np.arange(1,n+1)
    r=np.minimum.accumulate(r[::-1])[::-1]; q=np.empty(n); q[o]=np.minimum(r,1); return q
