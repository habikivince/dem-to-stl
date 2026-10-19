import subprocess, sys, resource, time, os
script, out, *extra = sys.argv[1:]
cmd=['/home/claude/venv/bin/python', script,'--input','/home/claude/bench/dem.tif','--cx','503000','--cy','5003000','--radius','2900','--diameter-mm','250','--vexag','1.5','--base-mm','3','--print-spacing-mm','0.2','--output',out,*extra]
t=time.time(); r=subprocess.run(cmd,capture_output=True,text=True); dt=time.time()-t
ru=resource.getrusage(resource.RUSAGE_CHILDREN)
print(f"{dt:.2f} {ru.ru_maxrss/1e3:.0f} rc={r.returncode}")
