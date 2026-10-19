import numpy as np, glob, os, sys
from osgeo import gdal
gdal.UseExceptions()
FILES=sorted(glob.glob('/home/claude/bench/real/*.tif*'))
def load(p):
    ds=gdal.Open(p); b=ds.GetRasterBand(1); a=b.ReadAsArray().astype(np.float64); nd=b.GetNoDataValue()
    bad=~np.isfinite(a)
    if nd is not None: bad|=(a<=nd+1e-6*abs(nd)) if nd<-1e6 else (a==nd)
    a[bad]=np.nan; return a
def surf(za,zb,zc,zd,u,v,ac):
    t_ac=np.where(v<=u, za+u*(zb-za)+v*(zc-zb), za+v*(zd-za)+u*(zc-zd))
    t_bd=np.where(u+v<=1, za+u*(zb-za)+v*(zd-za), zc+(1-u)*(zd-zc)+(1-v)*(zb-zc))
    return np.where(ac,t_ac,t_bd)
def run(z,s):
    n=z.shape[0]; m=(n-1)//s
    A=z[0:m*s:s,0:m*s:s]; B=z[0:m*s:s,s:m*s+1:s]; C=z[s:m*s+1:s,s:m*s+1:s]; D=z[s:m*s+1:s,0:m*s:s]
    rules={'fixe a-c':np.ones_like(A,bool),'fixe b-d':np.zeros_like(A,bool),'adaptative':np.abs(A-C)<=np.abs(B-D)}
    err={k:0.0 for k in rules}; err['oracle']=0.0; cnt=0
    for dv in range(s+1):
        for du in range(s+1):
            u,v=du/s,dv/s; truth=z[dv:dv+m*s:s,du:du+m*s:s][:m,:m]
            e1=(surf(A,B,C,D,u,v,rules['fixe a-c'])-truth)**2; e2=(surf(A,B,C,D,u,v,rules['fixe b-d'])-truth)**2; e3=(surf(A,B,C,D,u,v,rules['adaptative'])-truth)**2
            err['fixe a-c']+=e1.sum(); err['fixe b-d']+=e2.sum(); err['adaptative']+=e3.sum(); err['oracle']+=np.minimum(e1,e2).sum(); cnt+=m*m
    return {k:np.sqrt(e/cnt) for k,e in err.items()}
def windows(a,size,maxn=6):
    H,W=a.shape; out=[]
    for r in range(0,H-size+1,size//2):
        for c in range(0,W-size+1,size//2):
            w=a[r:r+size,c:c+size]
            if np.isfinite(w).all() and np.ptp(w)>5*np.finfo(float).eps and np.ptp(w)>0: out.append((np.ptp(w),r,c))
    out.sort(reverse=True)    # fenêtres les plus accidentées d'abord, sans recouvrement fort
    chosen=[]
    for rg,r,c in out:
        if all(abs(r-r2)>=size//2 or abs(c-c2)>=size//2 for _,r2,c2 in chosen): chosen.append((rg,r,c))
        if len(chosen)>=maxn: break
    return [(r,c,a[r:r+size,c:c+size]) for _,r,c in chosen]
if __name__=='__main__':
    rows=[]
    for f in FILES:
        a=load(f); name=os.path.basename(f)[:30]
        for s in (4,8):
            size=8*16*2+1 if s==8 else 4*16*4+1      # 257 pour s=8 et pour s=4 (257 = 4*64+1 = 8*32+1)
            size=257
            ws=windows(a,size)
            if not ws: continue
            agg={k:[] for k in ('fixe a-c','fixe b-d','adaptative','oracle')}
            for r,c,w in ws:
                res=run(w,s); base=res['fixe b-d']
                for k in agg: agg[k].append(res[k]/base)
            rows.append((name,s,len(ws),{k:np.mean(v) for k,v in agg.items()},{k:np.max(v) for k,v in agg.items()}))
            print(f"{name:30} pas={s} fenêtres={len(ws)}  moyenne / votre diagonale actuelle (b-d) : a-c {rows[-1][3]['fixe a-c']:.3f} | adaptative {rows[-1][3]['adaptative']:.3f} | oracle {rows[-1][3]['oracle']:.3f}   (pire fenêtre : adaptative {rows[-1][4]['adaptative']:.3f})",flush=True)
    import pickle; pickle.dump(rows,open('/home/claude/bench/diag_real.pkl','wb'))
