# Input files — docking, molecular dynamics, MM-GBSA

---

## `docking/` — AutoDock Vina

| File | Description |
|---|---|
| `receptor_ZEB2.pdbqt` | receptor (ZEB2 CZF construct) |
| `ligand_kEuT.pdbqt` | ZINCAG000000kEuT |
| `ligand_kH6D.pdbqt` | ZINCCJ000000kH6D |
| `ligand_BLM2.pdbqt` | ZINCxN000007BLM2 |
| `ligand_awc6.pdbqt` | ZINCzS000000awc6 |
| `ligand_aY0Z.pdbqt` | ZINCyD000000aY0Z |
| `vina_config.txt` | search box and sampling parameters |

`vina_config.txt`:
```
center_x = -5.744    center_y = 1.926    center_z = 1.763
size_x = 30.0        size_y = 30.0       size_z = 30.0
exhaustiveness = 32  num_modes = 10      energy_range = 3
seed = 42
```

---

## `gromacs/` — GROMACS

| File | Description |
|---|---|
| `topol.top`, `complex_topol.top` | system topology |
| `topol_Protein_chain_A.itp` | protein topology |
| `topol_Ion_chain_A2.itp`, `topol_Ion_chain_A3.itp` | two Zn²⁺ ions |
| `LIG.itp`, `lig_atomtypes.itp` | ligand parameters |
| `posre_*.itp` | positional restraints |
| `index.ndx` | index groups |

Force field: **amber14sb.ff**, water **TIP3P**.

---

## `gmx_mmpbsa/` — gmx_MMPBSA

`mmpbsa.in`:
```
&general
  startframe=1, endframe=500, interval=1,
  temperature=310.00, PBRadii=3,
/
&gb
  igb=8, saltcon=0.150,
/
&decomp
  idecomp=2, dec_verbose=3, print_res="within 4",
/
```
