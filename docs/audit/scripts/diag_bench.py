import numpy as np
def fbm(n,beta,seed,ridged=False):
    rng=np.random.default_rng(seed)
    f=np.fft.fftfreq(n)[:,None]**2+np.fft.fftfreq(n)[None,:]**2
    f[0,0]=1; amp=f**(-beta/4); amp[0,0]=0
    z=np.fft.ifft2(amp*np.fft.fft2(rng.standard_normal((n,n)))).real
    z=(z-z.mean())/z.std()
    if ridged: z=-np.abs(z)    # crêtes/vallées vives
    return z
def surf(za,zb,zc,zd,u,v,diag_ac):
    t_ac=np.where(v<=u, za+u*(zb-za)+v*(zc-zb), za+v*(zd-za)+u*(zc-zd))
    t_bd=np.where(u+v<=1, za+u*(zb-za)+v*(zd-za), zc+(1-u)*(zd-zc)+(1-v)*(zb-zc))
    return np.where(diag_ac,t_ac,t_bd)
def run(z,s):
    n=z.shape[0]; m=(n-1)//s
    A=z[0:m*s:s,0:m*s:s]; B=z[0:m*s:s,s:m*s+1:s]; C=z[s:m*s+1:s,s:m*s+1:s]; D=z[s:m*s+1:s,0:m*s:s]  # a=(r,c) b=(r,c+1) c=(r+1,c+1) d=(r+1,c)
    # coordonnées: u = colonne, v = ligne
    rule={'fixe a-c':np.ones_like(A,bool),'fixe b-d':np.zeros_like(A,bool),
          'adaptative (|za-zc|<=|zb-zd|)':np.abs(A-C)<=np.abs(B-D),
          'adaptative inverse (pire cas)':np.abs(A-C)>np.abs(B-D)}
    err={k:0.0 for k in rule}; cnt=0
    for dv in range(s+1):
        for du in range(s+1):
            u=du/s; v=dv/s
            truth=z[dv:dv+m*s:s,du:du+m*s:s][:m,:m]
            # points tronqués: aligne les tailles
            for k,r in rule.items():
                est=surf(A,B,C,D,u,v,r)
                err[k]+=np.sum((est-truth[:m,:m])**2)
            cnt+=m*m
    return {k:np.sqrt(e/cnt) for k,e in err.items()}
n=2049
for name,beta,ridged in [('fBm lisse (β=2.4)',2.4,False),('fBm rugueux (β=1.8)',1.8,False),('crêtes/vallées vives (|fBm| β=2.2)',2.2,True)]:
    print(f"\n== {name}")
    z=fbm(n,beta,1,ridged)
    for s in (4,16):
        r=run(z,s); base=r['fixe a-c']
        print(f"  pas={s}: "+' | '.join(f"{k}: {v/base:.3f}" for k,v in r.items())+f"   (RMSE fixe = {base:.4f} σ)")
