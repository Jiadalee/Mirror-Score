# Prospective D-peptide design protocol (AMR targets)

This protocol turns a target protein structure into a ranked D-peptide binder
library using Mirror-Score. It was prepared for two antimicrobial-resistance
(AMR) targets: **LasR** (P. aeruginosa quorum-sensing receptor; PDB 6MVM/2UV0)
and **LecB** (P. aeruginosa fucose-binding lectin; PDB 1OXC/1OVS).

## Overview

```
target PDB -> mirror (x-axis) -> QC -> hotspot map -> RFdiffusion binder
backbones -> ProteinMPNN sequences -> mirror back to D-peptides
-> Boltz-2 mirror-space cofold -> Mirror-Score ranking -> top-K library
```

## Step 1 — Mirror the target

Reflect the target structure through the x-axis (Garton et al., PNAS 2018
convention) so that a D-peptide binder of the mirror-image target corresponds
to an L-peptide binder of the mirror-image protein, which structure-prediction
and design tools handle natively.

```python
from mirrorscore.mirror import mirror_pdb, mirror_qc
mirror_pdb("data/amr_targets/structures/LasR.pdb",
           "data/amr_targets/mirrored/LasR_Dform.pdb")
mirror_qc("data/amr_targets/mirrored/LasR_Dform.pdb")  # |phi_L+phi_D| ~ 0
```

Mirrored structures and QC metrics for both targets ship in
`data/amr_targets/` (LasR: residues 7-168, chain A; LecB: residues 1-114).

## Step 2 — Hotspot selection

Hotspots are target residues contacting the co-crystallized ligand:

- **LasR**: 27 residues from the 3-oxo-C12-HSL ligand (K4G), e.g.
  `A36,A38,A39,A40,A47,A50,A52,A56,A60,A61,A64,A70,A73,A75,A76,A79,A80,A88,
  A93,A101,A105,A110,A115,A125,A126,A127,A129`
- **LecB**: 11 residues from the fucose ligand (FUC):
  `A21,A22,A23,A45,A95,A96,A97,A99,A101,A103,A104`

Hotspot maps ship in `data/amr_targets/rfdiffusion_input/design_targets.json`.

## Step 3 — Binder backbone generation (RFdiffusion)

Generate L-peptide binder backbones against the mirrored target (GPU required;
one job per target, ~30 min for 50 designs):

```bash
python /app/RFdiffusion/scripts/run_inference.py \
    inference.input_pdb=/input/LasR_Dform.pdb \
    'contigmap.contigs=[A7-168/0 15-30]' \
    'ppi.hotspot_res=[A36,A38,A39,A40,A47,A50,A52,A56,A60,A61,A64,A70,A73,A75,A76,A79,A80,A88,A93,A101,A105,A110,A115,A125,A126,A127,A129]' \
    inference.output_prefix=/output/binder \
    inference.num_designs=50 \
    inference.model_directory_path=/app/RFdiffusion/models
```

For LecB use `contigmap.contigs=[A1-114/0 15-30]` and the LecB hotspot list.

## Step 4 — Sequence design and mirror back

Run ProteinMPNN on each backbone (complex context: target + binder chain),
then reverse-map the designed L-peptide sequences to their D-enantiomers.

## Step 5 — Mirror-Score ranking

Cofold each designed D-peptide against the **native** target with Boltz-2 and
rank by interface pLDDT (primary), breaking ties with interface PAE. This is
the feature that achieved within-family LOO Spearman rho = 0.90 (p = 0.006)
on the benchmark with perfect top-3 recovery of the tightest binders, where
ProteinMPNN NLL fails (rho = 0.18).

```python
from mirrorscore.boltz_features import extract_features
# after: boltz predict design.yaml (target + candidate D-peptide)
feats = extract_features("predictions/candidate_0/", "candidate_0")
rank_key = feats["boltz_complex_iplddt"]   # higher = better
```

## Validation status

- Ranking module validated on 10 benchmark complexes with measured affinities
  (see `data/benchmark/calibration_summary.json`).
- RFdiffusion steps for LasR/LecB are fully specified but were not executed
  in the benchmarking session (compute-quota limit); all inputs and commands
  are provided for reproduction.
