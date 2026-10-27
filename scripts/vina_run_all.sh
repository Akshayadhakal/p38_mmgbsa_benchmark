#!/bin/bash
cd ~/benchmark_project/docking
echo "name,seed,vina_top" > vina_scores.csv
for seed in 1 2 3; do
  for f in ligands/*.pdbqt; do
    n=$(basename "$f" .pdbqt)
    vina --receptor receptor.pdbqt --ligand "$f" \
      --center_x 22.54 --center_y 10.62 --center_z 31.13 \
      --size_x 28 --size_y 20 --size_z 20 \
      --exhaustiveness 32 --num_modes 9 --seed $seed --cpu 8 \
      --out out/${n}_s${seed}.pdbqt > logs/${n}_s${seed}.log 2>&1
    s=$(awk '$1==1 {print $2; exit}' logs/${n}_s${seed}.log)
    echo "$n,$seed,$s" >> vina_scores.csv
  done
done
echo DONE >> vina_scores.csv
