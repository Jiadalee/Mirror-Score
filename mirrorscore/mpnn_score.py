"""ProteinMPNN sequence-backbone compatibility scoring in mirror space.

Reuses the ProteinMPNN implementation vendored in Mirror-Peptidizer
(Apache-2.0, https://github.com/dahuilangda/Mirror-Peptidizer) and the original
ProteinMPNN (Dauparas et al. 2022). Scores the (mirrored) L-peptide sequence on
the complex backbone with the target held fixed.
"""
import os
import sys
import warnings
import numpy as np

warnings.filterwarnings("ignore")

# Path to the Mirror-Peptidizer checkout containing ProteinMPNN/
MPNN_REPO = os.environ.get(
    "MIRROR_PEPTIDIZER_REPO",
    "/workspace/mirror_peptidizer/Mirror-Peptidizer-main",
)

D_CODES = {"DAL","DAR","DSG","DAS","DCY","DGN","DGL","DHI","DIL","DLE","DLY",
           "MED","DPN","DPR","DSN","DTH","DTR","DTY","DVA","DSP"}


def prep_complex_for_mpnn(src_pdb, dst_pdb, peptide_chain=None, mirror=True):
    """Strip waters/ligands/H, keep polymer records, optionally mirror, and
    promote D-residue HETATM records to ATOM after mirroring (they become
    standard L-names in mirror space)."""
    from mirrorscore.mirror import mirror_pdb_lines, L_TO_D, D_TO_L
    with open(src_pdb) as f:
        lines = f.readlines()
    if mirror:
        lines = mirror_pdb_lines(lines)
    out = []
    for line in lines:
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        resname = line[17:20].strip()
        resseq = line[22:26]
        # drop waters and non-polymer artifacts
        if resname in ("HOH", "WAT", "GOL", "SO4", "PO4", "CL", "NA", "EDO", "MPD", "ACT"):
            continue
        new_line = line
        if rec == "HETATM":
            # after mirroring, D-codes became L-codes; promote polymer residues to ATOM
            if resname in ("DAL","DAR","DSG","DAS","DCY","DGN","DGL","DHI","DIL","DLE","DLY",
                           "MED","DPN","DPR","DSN","DTH","DTR","DTY","DVA","DSP") or \
               resname in ("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
                           "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"):
                new_line = "ATOM  " + line[6:]
        out.append(new_line)
    with open(dst_pdb, "w") as f:
        f.writelines(out)
    return dst_pdb


def find_peptide_chain_original(pdb_path):
    """Identify the D-peptide chain in the ORIGINAL (unmirrored) file: the chain
    with the most D-residues (ties -> most residues). Chain IDs are preserved by
    mirroring, so the same ID is the L-peptide in mirror space."""
    from mirrorscore.features import load_chain_residues
    chains = load_chain_residues(pdb_path)
    best, best_key = None, (-1, -1)
    for c, reslist in chains.items():
        n_d = sum(1 for r in reslist if r.resname.strip().upper() in D_CODES)
        if n_d > 0:
            key = (n_d, len(reslist))
            if key > best_key:
                best, best_key = c, key
    if best is None:
        raise ValueError(f"no D-residue chain found in {pdb_path}")
    return best, best_key[1]


_MODEL = None
_UTILS = None


def _load():
    global _MODEL, _UTILS
    if _MODEL is not None:
        return _MODEL, _UTILS
    if MPNN_REPO not in sys.path:
        sys.path.insert(0, MPNN_REPO)
    os.chdir(MPNN_REPO)  # ProteinMPNN utils load weights relative to repo
    from utils import protein_mpnn as pm
    _UTILS = pm
    _MODEL = pm.load_model()
    return _MODEL, _UTILS


def mpnn_score(pdb_path, design_chain, score_mode="designed"):
    """ProteinMPNN NLL (lower = better backbone-sequence compatibility)."""
    model, pm = _load()
    import torch
    X, S, mask, chain_M, chain_M_pos, residue_idx, chain_encoding_all, _, _, _, _, _, _, _, _, _, _, _, _ = pm.prepare_inputs(pdb_path, design_chain=design_chain)
    with torch.no_grad():
        scores = pm.compute_native_score(model, X, S, mask, chain_M, chain_M_pos,
                                         residue_idx, chain_encoding_all, score_mode=score_mode)
    return float(np.asarray(scores).reshape(-1)[0])
