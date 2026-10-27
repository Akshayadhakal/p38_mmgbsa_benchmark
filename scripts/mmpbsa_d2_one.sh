#!/bin/bash
set -eo pipefail
D=${1:?need run dir}
trap 'echo "D2 FAILED: $D (line $LINENO)"' ERR
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mmpbsa
export PATH="$PATH:$HOME/miniconda3/envs/md/bin"
export GMXDATA="$HOME/miniconda3/envs/md/share/gromacs"
cd "$D"
if [ -f mmpbsa_d2.ok ]; then echo "$D: d2 already done"; exit 0; fi
for f in mm_C.xtc prod_sub.tpr topol_sub.top index_sub.ndx mmpbsa.ok; do
  if [ ! -f "$f" ]; then echo "$D: missing $f"; exit 1; fi
done
gnum() { awk -v g="$1" '/^\[/{if($2==g){print i; exit} i++}' index_sub.ndx; }
gsize() { awk -v g="$1" '/^\[/{s=($2==g)} !/^\[/ && s{n+=NF} END{print n+0}' index_sub.ndx; }
gp=$(gnum Protein); gl=$(gnum Other)
np=$(gsize Protein); nl=$(gsize Other); ns=$(gsize System)
if [ "$np" -ne 5637 ] || [ $((np+nl)) -ne "$ns" ]; then
  echo "$D: group sizes wrong (protein $np, ligand $nl, system $ns)"; exit 1
fi
cat > mmpbsa_d2.in << 'EOF'
&general
sys_name="run_d2",
startframe=251,
endframe=1001,
interval=3,
verbose=2,
/
&gb
igb=5,
saltcon=0.150,
intdiel=2.0,
/
EOF
mpirun -np 4 gmx_MMPBSA -O -i mmpbsa_d2.in -cs prod_sub.tpr -ci index_sub.ndx -cg $gp $gl \
  -ct mm_C.xtc -cp topol_sub.top -o final_results_d2.dat -eo final_results_d2.csv -nogui > mmpbsa_d2.out 2>&1
grep -q "^ΔTOTAL" final_results_d2.dat
touch mmpbsa_d2.ok
echo "$D d2 done: $(awk '/^Delta \(Complex/{f=1} f && /^ΔTOTAL/{print "dTOTAL", $2, "block SEM", $8; exit}' final_results_d2.dat)"
