import sys; sys.path.insert(0,'/home/claude/bench'); sys.path.insert(0,'/home/claude/impl/dem-to-stl')
import numpy as np, tin_real, adaptive_top
import matplotlib; matplotlib.use('Agg')
from matplotlib.tri import Triangulation, LinearTriInterpolator
def deviation(z, quad, tol):
    H1,W1=quad.shape; X=np.arange(W1+1,dtype=float); Y=-np.arange(H1+1,dtype=float)
    tris=adaptive_top.adaptive_top_triangles(quad,X,Y,z,tol)
    pts=tris.reshape(-1,3); key=np.round(pts[:,:2]).astype(np.int64)
    uniq,inv=np.unique(key[:,0]*100000+key[:,1],return_inverse=True)
    vx=np.zeros(len(uniq)); vy=np.zeros(len(uniq)); vz=np.zeros(len(uniq))
    vx[inv]=pts[:,0]; vy[inv]=pts[:,1]; vz[inv]=pts[:,2]
    tri=Triangulation(vx,vy,inv.reshape(-1,3)); interp=LinearTriInterpolator(tri,vz)
    # nœuds appartenant à au moins une cellule du masque
    cp=np.pad(quad,1); inn=cp[0:-1,0:-1]|cp[0:-1,1:]|cp[1:,0:-1]|cp[1:,1:]
    ii,jj=np.nonzero(inn[:H1+1,:W1+1])
    val=interp(X[jj],Y[ii]); ok=~val.mask
    dev=np.abs(np.asarray(val[ok])-z[ii[ok],jj[ok]])
    return len(tris), dev.max(), np.percentile(dev,99.9), (~ok).sum()
if __name__=='__main__':
    files={'USACE 5 m':'USACE-WRIS_example.tif','SheepMtn 62 m':'ChHarding__SheepMtn.tif','Cotopaxi 31 m':'fgebhart_m_cotopaxi.tiff','Kauai 249 m':'hicatlab_D_kauai.tif'}
    print("Écart vertical maximal entre la surface adaptative et les altitudes des nœuds du DEM, en multiples de la tolérance")
    print(f"{'relevé':16}{'tol (mm)':>9}{'triangles':>11}{'écart max / tol':>17}{'p99,9 % / tol':>15}{'nœuds non évalués':>19}")
    for label,fn in files.items():
        d=tin_real.prepare('/home/claude/bench/real/'+fn,500)
        for tol in (0.02,0.05,0.1):
            n,mx,p,miss=deviation(d['z'],d['quad'],tol)
            print(f"{label:16}{tol:>9}{n:>11}{mx/tol:>17.2f}{p/tol:>15.2f}{miss:>19}",flush=True)
