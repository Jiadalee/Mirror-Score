# Mirror-Score: methods detail

## 1. Mirror-space convention

All mirror operations reflect atomic x-coordinates (x -> -x), following the
Garton et al. (PNAS 2018) convention used by Mirror-Peptidizer. A D-peptide /
L-protein complex is converted to an all-L complex by reflecting the entire
structure and relabeling D-residues (CCD codes DAL, DAR, DSG, ...) to their L
parents. Chirality QC verifies backbone torsion inversion: |phi_L + phi_D| ~ 0
and |psi_L + psi_D| ~ 0 for every residue.

## 2. Scoring features

For a D-peptide/L-protein crystal (or predicted) complex:

### 2.1 ProteinMPNN mirror-space NLL
The mirrored complex (all-L) is scored with ProteinMPNN (vendored from
Mirror-Peptidizer; Dauparas et al. 2022 weights v_48_020) in
`score_mode='designed'`: the negative log-likelihood of the native (mirrored)
peptide sequence given the complex backbone, with the target chain fixed.
Lower NLL = better backbone-sequence compatibility.

### 2.2 Interface descriptors
Computed from the original (non-mirrored) structure at 5.0 A heavy-atom
cutoff: atom-atom contact count, peptide/target residues in contact,
hydrogen bonds, salt bridges, peptide and target net charges, charge
complementarity, peptide hydrophobic fraction, peptide length.

### 2.3 Boltz-2 mirror-space confidence
The complex sequences (both chains as L-form) are submitted to Boltz-2
co-folding. Mirror-space confidence metrics are chirality-agnostic affinity
features: ipTM, min cross-chain ipTM, interface pLDDT, interface PAE
(mean PAE between peptide and target tokens), peptide self-PAE, complex PDE.

## 3. Benchmark

31 heterochiral D-peptide/L-protein crystal complexes from the PDB,
identified by systematic full-text search and verified by D-residue CCD codes
("D-PEPTIDE LINKING" type). Families: viral entry/gp41 (11), cancer PPI
(MDM2/CHIP, 13), angiogenesis (3), AMR (1), designed PPI (1), enzyme
substrate (1), antibody (1).

18 complexes have literature-verified affinities (KD or IC50) and form the
calibration table. Affinity provenance is recorded per entry in
`data/benchmark/benchmark_complexes.csv`.

## 4. Calibration protocol

- **Within-family LOO**: ridge regression (standardized features) from
  features to log10 KD, leave-one-out within each family. This matches the
  realistic deployment scenario (ranking a design series against one target).
- **Leave-one-family-out (LOFO)**: tests cross-family transfer. Reported
  honestly: with current features, cross-family transfer fails (sign flips
  between families), which is itself a key finding.
- **L-transfer check**: family-matched L-peptide complexes (PMI-MDM2,
  17-28p53-MDM2, pDIQ-MDM2, Hsc70-EEVD-CHIP, TE33 epitope) verify that
  feature-affinity relationships are chirality-invariant.

## 5. Prospective AMR design workflow

1. Mirror D-form targets: LasR (PDB 6MVM), LecB (PDB 1OXC).
2. Hotspot selection from co-crystallized ligands (LasR: K4G autoinducer;
   LecB: fucose) - any target residue within 4.5 A of the ligand.
3. RFdiffusion binder design on the mirrored target
   (`contigmap.contigs=[A7-168/0 15-30]`, hotspot-constrained).
4. ProteinMPNN sequence design on backbones.
5. Mirror sequences back to D-form; rank with Mirror-Score.
