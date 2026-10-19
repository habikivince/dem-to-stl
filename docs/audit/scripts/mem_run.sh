#!/bin/bash
# usage: mem_run.sh sortie [options...]   — plafonne l'espace d'adressage à 3,6 Go pour ne pas faire tuer le conteneur
out=$1; shift
ulimit -v 3600000
/home/claude/venv/bin/python /home/claude/bench/time_run.py /home/claude/impl/dem-to-stl/circle_dem_to_stl.py "$out" --force "$@" 2>&1 | tail -3
