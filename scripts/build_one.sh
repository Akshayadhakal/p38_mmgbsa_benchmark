#!/bin/bash
set -euo pipefail
L=$1
ROOT=$HOME/benchmark_project
REF=$ROOT/md/lig_p38a_3fly
W=$ROOT/md/$L
trap 'echo "BUILD FAILED: $L (line $LINENO)"' ERR

mkdir -p "$W"; cd "$W"
if [ -f build.ok ]; then echo "$L: already built"; exit 0; fi

# 1. ligand parameters: GAFF2, AM1-BCC, net charge 0
cp "$ROOT/docking/ligands/$L.mol" lig.mol
obabel lig.mol -O lig.mol2 > obabel.log 2>&1
python $HOME/benchmark_project/scripts/fix_mol2.py lig.mol2
rm -rf LIG.acpype
acpype -i lig.mol2 -b LIG -c bcc -n 0 -a gaff2 > acpype.out 2>&1
test -s LIG.acpype/LIG_GMX.itp
q=$(awk '/\[ atoms \]/{f=1;next} /\[ bonds \]/{f=0} f && $1 ~ /^[0-9]+$/ {s+=$7} END{printf "%.3f", s}' LIG.acpype/LIG_GMX.itp)
awk -v q="$q" 'BEGIN{exit (q>0.01 || q<-0.01)}'

# ligand must sit in the same pocket as the protein file
python - << 'EOF'
import math, sys
pts = [(float(l[20:28]), float(l[28:36]), float(l[36:44])) for l in open('LIG.acpype/LIG_GMX.gro').read().splitlines()[2:-1]]
c = [sum(p[i] for p in pts) / len(pts) for i in range(3)]
d = math.dist(c, (2.254, 1.062, 3.113))
print("ligand center distance to pocket: %.2f nm" % d)
sys.exit(0 if d < 1.0 else 1)
EOF

# 2. protein topology, identical for every ligand (histidine states from the benchmark file)
cp "$REF/prot_fix.pdb" .
printf "1\n1\n1\n1\n1\n1\n1\n2\n0\n0\n1\n0\n" | gmx pdb2gmx -f prot_fix.pdb -o prot.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh -his > pdb2gmx.log 2>&1
grep -q "Total charge -7.000" pdb2gmx.log

# 3. combine protein and ligand
cp LIG.acpype/LIG_GMX.itp LIG.itp
cp LIG.acpype/posre_LIG.itp .
printf '\n#ifdef POSRES\n#include "posre_LIG.itp"\n#endif\n' >> LIG.itp
name=$(awk '/\[ *moleculetype *\]/{f=1;next} f && NF && $1!~/^;/{print $1;exit}' LIG.itp)
[ "$name" = "LIG" ]
n1=$(sed -n 2p prot.gro | tr -d ' ')
n2=$(sed -n 2p LIG.acpype/LIG_GMX.gro | tr -d ' ')
box=$(tail -1 prot.gro)
{ echo complex; echo $((n1+n2)); sed '1,2d;$d' prot.gro; sed '1,2d;$d' LIG.acpype/LIG_GMX.gro; echo "$box"; } > complex.gro
sed -i '/forcefield.itp"/a #include "LIG.itp"' topol.top
printf "%-20s 1\n" "$name" >> topol.top

# 4. box, water, ions
gmx editconf -f complex.gro -o box.gro -c -d 1.0 -bt dodecahedron > editconf.log 2>&1
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top > solvate.log 2>&1
cat > ions.mdp << 'EOF'
integrator = steep
nsteps = 0
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
EOF
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1 > grompp_ions.log 2>&1
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 > genion.log 2>&1
na=$(awk '/\[ molecules \]/{f=1} f && $1=="NA"{print $2}' topol.top)
cl=$(awk '/\[ molecules \]/{f=1} f && $1=="CL"{print $2}' topol.top)
[ $((na-cl)) -eq 7 ]

# 5. minimize (no warnings allowed)
cp "$REF/em.mdp" .
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr > grompp_em.log 2>&1
gmx mdrun -deffnm em -nb gpu > em.out 2>&1
grep -q "converged to Fmax" em.log

echo "$L built: protein $n1 + ligand $n2 atoms, NA $na CL $cl, ligand charge $q"
grep "Maximum force" em.log
touch build.ok
