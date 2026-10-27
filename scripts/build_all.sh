#!/bin/bash
cd ~/benchmark_project/md
for f in ~/benchmark_project/docking/ligands/*.mol; do
  L=$(basename "$f" .mol)
  [ "$L" = "lig_p38a_3fly" ] && continue
  bash ~/benchmark_project/scripts/build_one.sh "$L" || echo "$L" >> build_failed.txt
done
echo BUILD_DONE
