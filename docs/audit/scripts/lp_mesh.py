import sys, os, io, contextlib; sys.path.insert(0,'/home/claude/impl/dem-to-stl')
from line_profiler import LineProfiler
import mesh, circle_dem_to_stl as c
lp=LineProfiler(); lp.add_function(mesh.build_mesh_from_arrays)
f=lp(c.build_mesh)
dev=open(os.devnull,'w')
with contextlib.redirect_stdout(dev), contextlib.redirect_stderr(dev):
    f('/home/claude/bench/dem.tif', None, None, 503000.0, 5003000.0, 2900.0, 250.0, 1.5, 3.0, 0.2, bottom=sys.argv[1])
st=lp.get_stats(); unit=st.unit
src=open('/home/claude/impl/dem-to-stl/mesh.py',encoding='utf-8').read().split('\n')
tot=0; rows=[]
for (fn,ln,name),timings in st.timings.items():
    if name!='build_mesh_from_arrays': continue
    for lineno,hits,t in timings:
        rows.append((t*unit,lineno,hits)); tot+=t*unit
print(f"build_mesh_from_arrays [{sys.argv[1]}] : {tot:.2f} s au total ; lignes > 0,12 s :")
for t,lineno,hits in sorted(rows,reverse=True)[:14]:
    if t>0.12: print(f"  {t:5.2f} s  l.{lineno:<4} ({hits} appels)  {src[lineno-1].strip()[:95]}")
