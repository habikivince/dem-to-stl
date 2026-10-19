import sys, os, subprocess, time, resource
sys.path.insert(0,'/home/claude/bench'); sys.path.insert(0,'/home/claude/impl/dem-to-stl')
import tin_real
p='/home/claude/bench/real/USACE-WRIS_example.tif'; d=tin_real.prepare(p,500)
common=['--input',p,'--lat',repr(d['lat']),'--lon',repr(d['lon']),'--radius',repr(d['R']),'--diameter-mm','250','--vexag',repr(d['vexag']*1.25),'--base-mm','3','--print-spacing-mm',sys.argv[1],'--bottom','merged','--force']
print(f"USACE 5 m, Ø250 mm, pas {sys.argv[1]} mm, relief ≈ 25 mm, fond merged")
for name,extra in (('dessus grille',[]),('adaptatif 0,05 mm',['--top','adaptive','--top-tol-mm','0.05']),('adaptatif 0,1 mm',['--top','adaptive','--top-tol-mm','0.1'])):
    o=f'/home/claude/bench/big_{name.split()[0]}_{name[-6:-3].replace(",","")}.stl'
    cmd=['bash','-c','ulimit -v 3600000; exec /home/claude/venv/bin/python /home/claude/bench/time_run.py /home/claude/impl/dem-to-stl/circle_dem_to_stl.py "$0" "$@"',o]
    # time_run.py fixe le DEM synthétique : on appelle donc le script directement et on mesure ici
    t=time.time(); r=subprocess.run(['bash','-c','ulimit -v 3600000; exec /home/claude/venv/bin/python /home/claude/impl/dem-to-stl/circle_dem_to_stl.py "$@"','x']+common+['--output',o]+extra,capture_output=True,text=True)
    dt=time.time()-t; ru=resource.getrusage(resource.RUSAGE_CHILDREN)
    if r.returncode: print(name,'ECHEC',r.stderr[-300:]); continue
    print(f"  {name:20} {dt:5.1f} s   {os.path.getsize(o)/1e6:7.1f} Mo   {(os.path.getsize(o)-84)/50:>10.0f} triangles",flush=True)
