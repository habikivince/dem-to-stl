import numpy as np
from scipy.ndimage import gaussian_filter
def planar_blocks(z, tol, smin=2, smax=64):
    """Quadtree sur la grille de sommets (N=2^k+1). Un bloc de s×s cellules est fusionné si
    l'écart max au plan des moindres carrés (s+1)² sommets <= tol. Renvoie (nb_blocs, nb_cellules_dans_blocs_fusionnés>1 cellule)."""
    N=z.shape[0]; n=N-1; blocks=0; merged_cells=0
    stack=[(r,c,smax) for r in range(0,n,smax) for c in range(0,n,smax)]
    s_idx={}
    def fit_ok(r,c,s):
        sub=z[r:r+s+1,c:c+s+1]
        yy,xx=np.mgrid[0:s+1,0:s+1]
        A=np.c_[xx.ravel(),yy.ravel(),np.ones((s+1)**2)]
        coef,res,_,_=np.linalg.lstsq(A,sub.ravel(),rcond=None)
        return np.abs(A@coef-sub.ravel()).max()<=tol
    while stack:
        r,c,s=stack.pop()
        if s<=smin or not fit_ok(r,c,s):
            if s<=smin:
                blocks+=1
                if s>1 and fit_ok(r,c,s): merged_cells+=s*s
                else: blocks+=s*s-1   # cellules unitaires isolées (s=2 non planaire -> 4 cellules)
            else:
                h=s//2; stack+= [(r,c,h),(r,c+h,h),(r+h,c,h),(r+h,c+h,h)]
        else:
            blocks+=1; merged_cells+=s*s
    return blocks,merged_cells
h=np.fromfile('realdem/Terrain2STL-HEAD/hgt_files/N44W070.hgt',dtype='>i2').reshape(1201,1201).astype(float)
full=512*512
print("Fusion planaire en quadtree (optimiste : 2 triangles par bloc, sans surcoût de raccord). Part des triangles de la grille complète.")
print("Tolérance en % de l'amplitude : 0,05 %   0,10 %   0,25 %   0,50 %")
for sg,label in [(1.0,'σ=1'),(2.0,'σ=2'),(4.0,'σ=4')]:
    for r0,c0 in [(0,0),(600,600)]:
        z=gaussian_filter(h,sg)[r0:r0+513,c0:c0+513]; R=np.ptp(z)
        out=[]
        for f in (0.0005,0.001,0.0025,0.005):
            b,_=planar_blocks(z,f*R); out.append(f"{100*b/full:6.1f} %")
        print(f"  SRTM lissé {label} fenêtre({r0},{c0}):  "+"  ".join(out))
