#!/bin/bash
cd /home/claude/bench
for f in "$@"; do s=$(date +%s)
  prusa-slicer --export-gcode --layer-height 0.2 --first-layer-height 0.2 --fill-density 10% --bed-shape 0x0,250x0,250x210,0x210 --output $f.gcode $f.stl > $f.slice.log 2>&1
  echo "rc=$? durée=$(( $(date +%s)-s ))s" >> $f.slice.log
done
