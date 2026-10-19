import sys, subprocess
sys.path.insert(0,'/home/claude/bench'); sys.path.insert(0,'/home/claude/impl/dem-to-stl')
import tin_real
p='/home/claude/bench/real/USACE-WRIS_example.tif'; d=tin_real.prepare(p,500)
common=['--input',p,'--lat',repr(d['lat']),'--lon',repr(d['lon']),'--radius',repr(d['R']),'--diameter-mm','70','--vexag',repr(d['vexag']*0.5),'--base-mm','2','--print-spacing-mm','0.2','--bottom','merged','--force']
for name,extra in (('s_grid',[]),('s_tin05',['--top','adaptive','--top-tol-mm','0.05']),('s_tin10',['--top','adaptive','--top-tol-mm','0.1'])):
    r=subprocess.run(['/home/claude/venv/bin/python','/home/claude/impl/dem-to-stl/circle_dem_to_stl.py']+common+['--output',f'/home/claude/bench/{name}.stl']+extra,capture_output=True,text=True)
    print(name,'rc',r.returncode)
