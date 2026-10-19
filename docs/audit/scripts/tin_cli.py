import sys, os, subprocess, time, resource
sys.path.insert(0,'/home/claude/bench'); sys.path.insert(0,'/home/claude/impl/dem-to-stl')
import tin_real, check_stl, numpy as np
files={'USACE 5 m':'USACE-WRIS_example.tif','SheepMtn 62 m':'ChHarding__SheepMtn.tif','Cotopaxi 31 m':'fgebhart_m_cotopaxi.tiff','Kauai 249 m':'hicatlab_D_kauai.tif'}
out=[]
for label,fn in files.items():
    p='/home/claude/bench/real/'+fn; d=tin_real.prepare(p,500)
    base=['/home/claude/venv/bin/python','/home/claude/impl/dem-to-stl/circle_dem_to_stl.py','--input',p,'--lat',repr(d['lat']),'--lon',repr(d['lon']),'--radius',repr(d['R']),
          '--diameter-mm','200','--vexag',repr(d['vexag']),'--base-mm','3','--print-spacing-mm','0.4','--bottom','merged']
    row=[label]
    for name,extra in (('grille',[]),('adaptatif 0,02',['--top','adaptive','--top-tol-mm','0.02']),('adaptatif 0,05',['--top','adaptive','--top-tol-mm','0.05']),('adaptatif 0,1',['--top','adaptive','--top-tol-mm','0.1'])):
        o=f'/tmp/tc_{abs(hash(label+name))%10**6}.stl'; t=time.time()
        r=subprocess.run(base+['--output',o]+extra,capture_output=True,text=True); dt=time.time()-t
        if r.returncode: print(label,name,'ECHEC',r.stdout[-200:],r.stderr[-300:]); continue
        tr=check_stl.read_stl(o); res={n:ok for n,ok,_ in check_stl.check(tr,200.0,3.0,tol_mm=2.0)}
        row.append((name,len(tr),all(res.values()),dt,os.path.getsize(o)/1e6,o))
    out.append(row)
    g=row[1][1]
    print(f"{label:14} relief 20 mm, pixel de travail {d['pix']:.1f} m")
    for name,n,ok,dt,mb,o in row[1:]: print(f"    {name:16} {n:>8} triangles ({100*n/g:5.1f} % de la grille)  {mb:6.1f} Mo  {dt:4.1f} s  contrôles check_stl: {'OK' if ok else 'ÉCHEC'}",flush=True)
import pickle; pickle.dump([(r[0],[(a,b,c,dd,e,f) for a,b,c,dd,e,f in r[1:]]) for r in out],open('/home/claude/bench/tin_cli.pkl','wb'))
