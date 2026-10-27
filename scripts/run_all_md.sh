#!/bin/bash
cd ~/benchmark_project/md
while [ ! -f lig_p38a_3fly/prod.gro ]; do sleep 60; done
for r in $(seq $1 $2); do
  for f in ~/benchmark_project/docking/ligands/*.mol; do
    L=$(basename "$f" .mol)
    if [ "$L" = "lig_p38a_3fly" ] && [ "$r" -eq 1 ]; then continue; fi
    bash ~/benchmark_project/scripts/run_rep.sh "$L" "$r" >> md_all.log 2>&1 || echo "$L rep$r" >> md_failed.txt
  done
done
echo ALLDONE >> md_all.log
