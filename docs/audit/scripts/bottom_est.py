import numpy as np, time
def cover(mask):
    ny,nx=mask.shape; todo=mask.copy(); rects=[]
    for j in range(ny):
        i=0
        while i<nx:
            if not todo[j,i]: i+=1; continue
            i1=i
            while i1+1<nx and todo[j,i1+1]: i1+=1
            j1=j
            while j1+1<ny and todo[j1+1,i:i1+1].all(): j1+=1
            todo[j:j1+1,i:i1+1]=False; rects.append((i,i1,j,j1)); i=i1+1
    return rects
def estimate(diam_mm,sp):
    n=int(round(diam_mm/sp)); yy,xx=np.mgrid[0:n,0:n]; c=(n-1)/2
    cells=((xx+0.5-n/2)**2+(yy+0.5-n/2)**2)<=(n/2)**2
    cells=cells[:-1,:-1] if False else cells
    rects=cover(cells)
    merged=np.zeros_like(cells)
    for i0,i1,j0,j1 in rects: merged[j0:j1+1,i0:i1+1]=True
    P=np.pad(merged,1)                            # nœud (j,i) = coin SO de la cellule (j,i): entouré par cellules (j-1..j, i-1..i)
    # req[j,i] sur grille de nœuds (ny+1,nx+1): vrai si les 4 cellules autour ne sont pas toutes fusionnées
    ny,nx=cells.shape
    Q=np.pad(merged,1,constant_values=False)
    req=~(Q[1:,1:]&Q[1:,:-1]&Q[:-1,:-1]&Q[:-1,1:])
    for i0,i1,j0,j1 in rects:
        req[j0,i0]=req[j0,i1+1]=req[j1+1,i1+1]=req[j1+1,i0]=True
    tri=0; worst=0
    for i0,i1,j0,j1 in rects:
        ring=req[j0,i0:i1+2].sum()+req[j0+1:j1+2,i1+1].sum()+req[j1+1,i0:i1+1].sum()+req[j0+1:j1+1,i0].sum()
        tri+=ring
    full=2*int(cells.sum())
    return n,int(cells.sum()),full,len(rects),int(tri)
for d,s in [(250,0.2),(120,0.2),(250,0.1)]:
    t=time.time(); n,c,full,nr,tri=estimate(d,s)
    print(f"Ø{d} mm, pas {s} mm: grille {n}², {c:,} cellules → fond actuel = {full:,} triangles | couverture gloutonne: {nr} rectangles → fond fusionné ≈ {tri:,} triangles ({100*(1-tri/full):.2f} % en moins)  [{time.time()-t:.1f}s]")
    top=full; print(f"    total (dessus+fond, hors parois) : {2*full:,} → {full+tri:,}  (≈ −{100*(1-(full+tri)/(2*full)):.1f} %)")
