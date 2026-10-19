import sys, os, glob, time
sys.path.insert(0,'/home/claude/impl/dem-to-stl'); sys.setrecursionlimit(20000)
import numpy as np
from osgeo import gdal, osr
import raster, mesh
gdal.UseExceptions()

def rtin_masked(z, cells, tols):
    """RTIN sur z (N,N), N=2^k+1 ; cells (N-1,N-1) booléen = cellules du solide.
    Nœuds du bord du masque : erreur forcée à +inf (le contour en escalier est conservé exactement) ;
    nœuds sans cellule voisine dans le masque : erreur nulle. Renvoie {tol: nb de triangles dans le masque}."""
    N=z.shape[0]; E=np.zeros((N,N))
    c=np.pad(cells,1,constant_values=False)
    around=np.stack([c[1:,1:][:N,:N] if False else c[:-1,:-1],c[:-1,1:],c[1:,:-1],c[1:,1:]])    # 4 cellules autour de chaque nœud (N,N)
    inside_any=around.any(axis=0); inside_all=around.all(axis=0)
    forced=inside_any & ~inside_all; outside=~inside_any
    def err(ax,ay,bx,by,cx,cy):
        if abs(ax-cx)+abs(ay-cy)<=1: return 0.0
        mx=(ax+bx)>>1; my=(ay+by)>>1
        if forced[my,mx]: e=np.inf
        elif outside[my,mx]: e=0.0
        else: e=abs(z[my,mx]-(z[ay,ax]+z[by,bx])/2)
        e=max(e,err(cx,cy,ax,ay,mx,my),err(bx,by,cx,cy,mx,my))
        if e>E[my,mx]: E[my,mx]=e
        return e
    last=N-1
    err(last,last,0,0,last,0); err(0,0,last,last,0,last)
    cell=np.zeros((N,N),bool); cell[:N-1,:N-1]=cells
    def cnt(ax,ay,bx,by,cx,cy,tol):
        mx=(ax+bx)>>1; my=(ay+by)>>1
        if abs(ax-cx)+abs(ay-cy)<=1 or E[my,mx]<=tol:
            gx=(ax+bx+cx)/3.0; gy=(ay+by+cy)/3.0
            return 1 if cell[int(gy),int(gx)] else 0
        return cnt(cx,cy,ax,ay,mx,my,tol)+cnt(bx,by,cx,cy,mx,my,tol)
    out={}
    for t in tols: out[t]=cnt(last,last,0,0,last,0,t)+cnt(0,0,last,last,0,last,t)
    return out

def prepare(path, n_cells=500, diam_mm=200.0, relief_mm=20.0):
    ds=gdal.Open(path); b=ds.GetRasterBand(1); a=b.ReadAsArray().astype(float); nd=b.GetNoDataValue()
    bad=~np.isfinite(a)
    if nd is not None: bad|=(a==nd)|(a<-1e4)
    a[bad]=np.nan; H,W=a.shape; gt=ds.GetGeoTransform()
    # centre : fenêtre (0,7 x la plus petite dimension) au plus fort écart-type, entièrement valide à 98 %
    w=int(0.7*min(H,W)); best=None
    for r in np.linspace(0,H-w,3).astype(int):
        for c in np.linspace(0,W-w,3).astype(int):
            sub=a[r:r+w,c:c+w]
            if np.isfinite(sub).mean()<0.6: continue
            sd=np.nanstd(sub)
            if best is None or sd>best[0]: best=(sd,r,c)
    if best is None: return None
    _,r,c=best; pr,pc=r+w/2,c+w/2
    x=gt[0]+pc*gt[1]+pr*gt[2]; y=gt[3]+pc*gt[4]+pr*gt[5]
    src=osr.SpatialReference(); src.ImportFromWkt(ds.GetProjection()); src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    wgs=osr.SpatialReference(); wgs.ImportFromEPSG(4326); wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    lon,lat,_=osr.CoordinateTransformation(src,wgs).TransformPoint(x,y)
    wkt,cx,cy=raster.resolve_target_crs_and_center(path,lat,lon,None,None)
    px_m=raster.source_pixel_size_m(path,wkt,cx,cy)
    R=0.5*w*px_m*0.95                                      # rayon : presque la fenêtre choisie
    scale=(diam_mm/2)/R; pix=2*R/n_cells
    elev,xs,ys,nodata=raster.warp_window(path,wkt,cx,cy,R,pix,progress_cb=lambda *a:1,fill=True)
    valid=(elev!=nodata) if nodata is not None else np.ones_like(elev,bool)
    inside=(xs[None,:]-cx)**2+(ys[:,None]-cy)**2<=R**2
    _,quad=mesh.clean_mask(inside&valid)
    zmin=elev[valid&inside].min(); zmax=elev[valid&inside].max()
    vexag=relief_mm/((zmax-zmin)*scale)
    z=np.where(valid,(elev-zmin)*scale*vexag,0.0)
    return dict(z=z,quad=quad,R=R,px_src=px_m,pix=pix,relief_m=zmax-zmin,vexag=vexag,scale=scale,lat=lat,lon=lon)

def to_pow2(z,quad):
    H1,W1=quad.shape; n=max(H1+1,W1+1); N=2
    while N+1<n: N*=2
    N+=1
    zp=np.pad(z,((0,N-z.shape[0]),(0,N-z.shape[1])),mode='edge'); cp=np.pad(quad,((0,N-1-H1),(0,N-1-W1)))
    return zp,cp

if __name__=='__main__':
    n_cells=int(sys.argv[1]) if len(sys.argv)>1 else 500
    tols=[0.02,0.05,0.1,0.2]
    files=[f for f in sorted(glob.glob('/home/claude/bench/real/*.tif*')) if True] if len(sys.argv)<3 else sys.argv[2:]
    print(f"Grille de {n_cells} cellules sur Ø200 mm (pas {200/n_cells:.2f} mm), relief normalisé à 20 mm. Part des triangles du dessus conservée par RTIN à bord forcé (borne haute d'un TIN) :")
    print(f"{'relevé':32}{'pixel src':>10}{'cellules':>10}  "+"  ".join(f"{t:>5} mm" for t in tols))
    for f in files:
        d=prepare(f,n_cells)
        if d is None: print(os.path.basename(f)[:30],"aucune fenêtre valide"); continue
        zp,cp=to_pow2(d['z'],d['quad']); t0=time.time()
        r=rtin_masked(zp,cp,tols); full=2*int(d['quad'].sum())
        print(f"{os.path.basename(f)[:30]:32}{d['px_src']:>8.1f} m{int(d['quad'].sum()):>10}  "+"  ".join(f"{100*r[t]/full:>6.1f}%" for t in tols)+f"   [{time.time()-t0:.0f}s]",flush=True)
