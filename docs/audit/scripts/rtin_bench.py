import numpy as np, sys
sys.setrecursionlimit(10000)
def rtin_counts(z, tols):
    N=z.shape[0]; E=np.zeros((N,N))
    def err(ax,ay,bx,by,cx,cy):
        if abs(ax-cx)+abs(ay-cy)<=1: return 0.0      # triangle de la grille de pixels : pas de milieu entier
        mx=(ax+bx)>>1; my=(ay+by)>>1
        e=abs(z[my,mx]-(z[ay,ax]+z[by,bx])/2)
        e=max(e,err(cx,cy,ax,ay,mx,my),err(bx,by,cx,cy,mx,my))
        if e>E[my,mx]: E[my,mx]=e
        return e
    last=N-1
    err(last,last,0,0,last,0); err(0,0,last,last,0,last)
    out={}
    def cnt(ax,ay,bx,by,cx,cy,tol):
        mx=(ax+bx)>>1; my=(ay+by)>>1
        if abs(ax-cx)+abs(ay-cy)<=1 or E[my,mx]<=tol: return 1
        return cnt(cx,cy,ax,ay,mx,my,tol)+cnt(bx,by,cx,cy,mx,my,tol)
    for t in tols:
        out[t]=cnt(last,last,0,0,last,0,t)+cnt(0,0,last,last,0,last,t)
    return out
def noise_fbm(n,beta,seed):
    rng=np.random.default_rng(seed)
    f=np.fft.fftfreq(n)[:,None]**2+np.fft.fftfreq(n)[None,:]**2; f[0,0]=1
    amp=f**(-beta/4); amp[0,0]=0
    return np.fft.ifft2(amp*np.fft.fft2(rng.standard_normal((n,n)))).real
h=np.fromfile('realdem/Terrain2STL-HEAD/hgt_files/N44W070.hgt',dtype='>i2').reshape(1201,1201).astype(float)
fracs=[0.0005,0.001,0.0025,0.005]   # 0,05 % ; 0,1 % (défaut hmterrain) ; 0,25 % ; 0,5 % de l'amplitude
def report(name,z):
    z=z[:513,:513]; rng=z.max()-z.min(); full=2*512*512
    r=rtin_counts(z,[f*rng for f in fracs])
    print(f"  {name:44} "+"  ".join(f"{f*100:.2f}%: {100*r[f*rng]/full:5.1f}%" for f in fracs))
print("Part des triangles de la grille complète conservée par RTIN (borne haute d'un TIN glouton), erreur max = x % de l'amplitude du relief")
print("SRTM3 Maine (réel, 513²)")
for r0,c0 in [(0,0),(0,600),(600,0),(600,600)]: report(f"  fenêtre ({r0},{c0}), relief {np.ptp(h[r0:r0+513,c0:c0+513]):.0f} m",h[r0:r0+513,c0:c0+513])
print("Synthétique 513²")
for beta in (3.0,2.2,1.6):
    z=noise_fbm(513,beta,3); z=(z-z.min())/(z.max()-z.min()); report(f"  fBm lisse→rugueux (β={beta})",z)
print("SRTM + bruit gaussien (simule LiDAR bruité), fenêtre (0,0)")
z0=h[0:513,0:513]; R=np.ptp(z0); rng=np.random.default_rng(5)
for s in (0.0,0.0005,0.001,0.002):
    report(f"  bruit σ = {s*100:.2f} % de l'amplitude",z0+rng.normal(0,s*R,z0.shape))
