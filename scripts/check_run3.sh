#!/bin/bash
# usage: check_run3.sh <run_dir> <index_file> <group_with_protein_and_ligand>
set -euo pipefail
trap 'echo "CHECK FAILED (line $LINENO)"' ERR
D=$1; NDX=$2; G=$3
S=$HOME/benchmark_project/scripts
cd "$D"
L=check3.log; : > $L
echo "$G" | gmx convert-tpr -s prod.tpr -n "$NDX" -o prod_sub.tpr >> $L 2>&1
printf "System\n" | gmx trjconv -s prod_sub.tpr -f prod.xtc -pbc whole -o mm_whole.xtc >> $L 2>&1
echo q | gmx make_ndx -f prod_sub.tpr -o index_sub.ndx >> $L 2>&1
printf "Protein\nSystem\n" | gmx trjconv -s prod_sub.tpr -f mm_whole.xtc -n index_sub.ndx -pbc mol -center -ur compact -o mm_B.xtc >> $L 2>&1
printf "System\n" | gmx trjconv -s prod_sub.tpr -f mm_B.xtc -dump 0 -o mm_B_f0.gro >> $L 2>&1
printf "System\n" | gmx trjconv -s mm_B_f0.gro -f mm_B.xtc -o mm_C.xtc -pbc nojump >> $L 2>&1
rm -f mm_whole.xtc mm_B.xtc
awk '/\[ molecules \]/{m=1} m && ($1=="SOL"||$1=="NA"||$1=="CL"){next} {print}' topol.top > topol_sub.top
printf "System\n" | gmx trjconv -s prod_sub.tpr -f mm_C.xtc -skip 50 -o mm_C_check.pdb >> $L 2>&1
scan=$(python3 $S/bondscan_multi.py topol_sub.top mm_C_check.pdb)
echo "bond scan: $scan"
echo "$scan" | grep -q "A: 0 |"
gmx distance -s prod_sub.tpr -f mm_C.xtc -n index_sub.ndx -nopbc \
  -select 'com of group "Protein" plus com of group "Other"' -oall com_C_nopbc.xvg >> $L 2>&1
dist=$(awk '!/^[#@]/{c++; if($2>3)b++; if(m==""||$2<m)m=$2; if($2>x)x=$2} END{printf "no-PBC COM distance (nm): min %.2f max %.2f | frames above 3 nm: %d of %d", m, x, b+0, c}' com_C_nopbc.xvg)
echo "$dist"
echo "$dist" | grep -q "above 3 nm: 0 of 1001"
printf "Backbone\nBackbone\n" | gmx rms -s mm_B_f0.gro -f mm_C.xtc -n index_sub.ndx -nomw -o rmsd_bb3.xvg -tu ns >> $L 2>&1
printf "Backbone\nOther\n" | gmx rms -s mm_B_f0.gro -f mm_C.xtc -n index_sub.ndx -nomw -o rmsd_lig3.xvg -tu ns >> $L 2>&1
printf "Protein\nOther\n" | gmx mindist -s prod_sub.tpr -f prod.xtc -od mindist3.xvg >> $L 2>&1
s() { awk '!/^[#@]/{if(c==0)f=$2; t+=$2; if($2>m)m=$2; c++; l=$2} END{printf "first %.3f mean %.3f max %.3f last %.3f (%d frames)", f, t/c, m, l, c}' "$1"; }
echo "backbone RMSD (nm): $(s rmsd_bb3.xvg)"
echo "ligand RMSD   (nm): $(s rmsd_lig3.xvg)"
echo "min distance  (nm): $(s mindist3.xvg)"
