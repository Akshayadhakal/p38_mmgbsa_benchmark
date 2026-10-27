# MM/GBSA protocol (fixed before any comparison with experiment)
- Frames: 5-20 ns of each 20 ns run, every 3rd frame (251 frames)
- gmx_MMPBSA 1.7.0, igb=5, saltcon=0.150, intdiel=1, no entropy
- Replicas: values averaged per ligand; each replica also analysed alone
- Sensitivity (secondary): intdiel=2
- Primary outcome: Spearman correlation with experimental dG, all 29 ligands, versus Vina (paired bootstrap)
- Secondary: without lig_p38a_2ee and lig_p38a_2l (1 nM); replica 1 and 2 separately; correlation with ligand heavy-atom count (exploratory)
Mon Oct 19 07:39:05 UTC 2026
