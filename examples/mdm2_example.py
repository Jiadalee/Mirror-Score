"""Worked example: rank the MDM2 D-peptide series (DPMI family) with
Mirror-Score features.

Requires the Mirror-Peptidizer checkout for vendored ProteinMPNN:
    export MIRROR_PEPTIDIZER_REPO=/path/to/Mirror-Peptidizer

Benchmark structures are in data/benchmark/structures/.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mirrorscore.mirror import mirror_pdb_file, mirror_qc
from mirrorscore.features import interface_descriptors
from mirrorscore.mpnn_score import mpnn_score, prep_complex_for_mpnn

# MDM2 D-peptide series with verified affinities
SERIES = [
    ("3TPX", "B", 0.22),   # DPMI-delta, 220 pM
    ("3IWY", "B", 53.0),   # (D)PMI-gamma, 53 nM
    ("3LNJ", "F", 219.0),  # DPMI-alpha, 219 nM
    ("8F0Z", None, 5500.0),  # H101, 5.5 uM
    ("8F10", None, 2400.0),  # H102, 2.4 uM
    ("8F12", None, 880.0),   # H103, 880 nM
]

STRUCT_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "benchmark", "structures")
PREP_DIR = os.path.join(os.path.dirname(__file__), "prep")
os.makedirs(PREP_DIR, exist_ok=True)


def main():
    from mirrorscore.mpnn_score import find_peptide_chain_original

    rows = []
    for pdb_id, chain_hint, kd in SERIES:
        src = os.path.join(STRUCT_DIR, f"{pdb_id}.pdb")
        if not os.path.exists(src):
            print(f"{pdb_id}: structure not found, skipping")
            continue

        # peptide chain: use hint or auto-detect (chain with most D-residues)
        chain = chain_hint or find_peptide_chain_original(src)

        # mirror-space MPNN NLL
        prep = os.path.join(PREP_DIR, f"{pdb_id}_mirror.pdb")
        prep_complex_for_mpnn(src, prep, peptide_chain=chain, mirror=True)
        nll = mpnn_score(prep, design_chain=chain, score_mode="designed")

        # interface descriptors (original chirality)
        desc = interface_descriptors(src, chain)

        rows.append((pdb_id, kd, nll, desc["n_atom_contacts"]))
        print(f"{pdb_id}: KD={kd:>8.1f} nM  NLL={nll:.3f}  contacts={desc['n_atom_contacts']}")

    # within-family ranking check (Spearman)
    from scipy.stats import spearmanr
    import numpy as np

    kds = np.log10([r[1] for r in rows])
    nlls = [r[2] for r in rows]
    rho, p = spearmanr(nlls, kds)
    print(f"\nSpearman(NLL, log10 KD) = {rho:.3f} (p={p:.3f})")
    print("Lower NLL should mean lower KD (tighter binding).")


if __name__ == "__main__":
    main()
