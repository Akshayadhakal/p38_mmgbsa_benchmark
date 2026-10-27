#!/bin/bash
set -euo pipefail
L=$1; R=$2
ROOT=$HOME/benchmark_project
REF=$ROOT/md/lig_p38a_3fly
W=$ROOT/md/$L
trap 'echo "RUN FAILED: $L rep$R (line $LINENO)"' ERR

cd "$W"
if [ "$L" != "lig_p38a_3fly" ]; then [ -f build.ok ]; fi
[ -f em.gro ]
mkdir -p rep$R; cd rep$R
if [ -f prod.gro ]; then echo "$L rep$R: already done"; exit 0; fi

n1=$(sed -n 2p ../prot.gro | tr -d ' ')
n2=$(sed -n 2p ../LIG.acpype/LIG_GMX.gro | tr -d ' ')
cp ../topol.top ../LIG.itp ../posre.itp ../posre_LIG.itp .
cp "$REF/nvt.mdp" "$REF/npt.mdp" .

# index with a "Prot_Lig" group (protein atoms, then ligand atoms)
if [ ! -f ../index_all.ndx ]; then
  echo q | gmx make_ndx -f ../em.gro -o ../index_all.ndx > make_ndx.log 2>&1
  awk -v N=$((n1+n2)) 'BEGIN{print ""; print "[ Prot_Lig ]"; for(i=1;i<=N;i++){printf "%d ", i; if(i%15==0) printf "\n"} printf "\n"}' >> ../index_all.ndx
fi
tot=$(awk '/^\[/{g=$2} g=="Prot_Lig" && !/^\[/ {n+=NF} END{print n}' ../index_all.ndx)
[ "$tot" -eq $((n1+n2)) ]
cp ../index_all.ndx .

NS=${NSTEPS:-10000000}
cat > prod.mdp << EOF
integrator              = md
dt                      = 0.002
nsteps                  = $NS
nstlog                  = 10000
nstenergy               = 5000
nstxout-compressed      = 10000
compressed-x-grps       = Prot_Lig
continuation            = yes
gen_vel                 = no
constraints             = h-bonds
constraint_algorithm    = lincs
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
fourierspacing          = 0.12
DispCorr                = EnerPres
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
EOF

gmx grompp -f nvt.mdp -c ../em.gro -r ../em.gro -p topol.top -o nvt.tpr > grompp_nvt.log 2>&1
gmx mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu > nvt.out 2>&1
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr > grompp_npt.log 2>&1
gmx mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu > npt.out 2>&1
gmx grompp -f prod.mdp -c npt.gro -t npt.cpt -p topol.top -n index_all.ndx -o prod.tpr > grompp_prod.log 2>&1
gmx mdrun -deffnm prod -nb gpu -pme gpu -bonded gpu > prod.out 2>&1
echo "$L rep$R done: $(grep -m1 Performance prod.log)"
