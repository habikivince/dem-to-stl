import subprocess, sys, resource, time, os
out, *extra = sys.argv[1:]
cmd=[ '/home/claude/venv/bin/python','/home/claude/impl/dem-to-stl/circle_dem_to_stl.py','--input','/home/claude/bench/dem.tif',
     '--cx','503000','--cy','5003000','--radius',os.environ.get('RADIUS','1400'),'--diameter-mm',os.environ.get('DIAM','120'),'--vexag','1.5','--base-mm','3',
     '--print-spacing-mm','0.2','--output',out,*extra]
t=time.time(); r=subprocess.run(cmd,capture_output=True,text=True)
ru=resource.getrusage(resource.RUSAGE_CHILDREN)
print(f"{os.path.basename(out)}: rc={r.returncode} durée={time.time()-t:.1f}s pic RAM={ru.ru_maxrss/1e3:.0f} Mo taille={os.path.getsize(out)/1e6:.1f} Mo")
if r.returncode: print(r.stdout[-600:], r.stderr[-600:])
