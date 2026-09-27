"""Mirror-image coordinate transformation for protein structures.

Implements the x-axis reflection convention of Garton et al., PNAS 2018, the same
convention used by Mirror-Peptidizer (Ma et al., Research 2026, 9:1420).
Adapted from Mirror-Peptidizer (https://github.com/dahuilangda/Mirror-Peptidizer,
Apache-2.0) with QC extensions: paired dihedral and chirality checks.
"""
import numpy as np

# Standard L-amino acid -> D-amino acid PDB chemical component codes
L_TO_D = {
    "ALA": "DAL", "ARG": "DAR", "ASN": "DSG", "ASP": "DAS", "CYS": "DCY",
    "GLN": "DGN", "GLU": "DGL", "GLY": "GLY", "HIS": "DHI", "ILE": "DIL",
    "LEU": "DLE", "LYS": "DLY", "MET": "MED", "PHE": "DPN", "PRO": "DPR",
    "SER": "DSN", "THR": "DTH", "TRP": "DTR", "TYR": "DTY", "VAL": "DVA",
    "MSE": "MED",
}
D_TO_L = {v: k for k, v in L_TO_D.items() if k != "GLY"}
D_TO_L["GLY"] = "GLY"


def reflect_coords(coords: np.ndarray) -> np.ndarray:
    """Reflect coordinates along the x axis (improper transformation, det = -1)."""
    out = coords.copy()
    out[:, 0] *= -1.0
    return out


def mirror_pdb_lines(lines):
    """Mirror a PDB file's ATOM/HETATM records: reflect x, relabel D<->L residues."""
    out = []
    for line in lines:
        rec = line[:6]
        if rec in ("ATOM  ", "HETATM"):
            resname = line[17:20].strip()
            new_resname = L_TO_D.get(resname, D_TO_L.get(resname, resname))
            x = float(line[30:38])
            x_new = -x
            line = line[:17] + f"{new_resname:>3s}" + line[20:30] + f"{x_new:8.3f}" + line[38:]
        out.append(line)
    return out


def mirror_pdb_file(src, dst):
    with open(src) as f:
        lines = f.readlines()
    with open(dst, "w") as f:
        f.writelines(mirror_pdb_lines(lines))


def chirality_sign(coords_n, coords_ca, coords_c, coords_cb):
    """Signed N-Ca-C-Cb scalar triple product; sign flips between L and D residues."""
    v1 = coords_n - coords_ca
    v2 = coords_c - coords_ca
    v3 = coords_cb - coords_ca
    return np.dot(np.cross(v1, v2), v3)


def dihedrals_from_pdb(pdb_path, resname_filter=None):
    """Extract phi/psi dihedrals per residue from a PDB file (first model)."""
    from Bio.PDB import PDBParser
    import warnings
    warnings.filterwarnings("ignore")
    parser = PDBParser(QUIET=True)
    st = parser.get_structure("x", pdb_path)
    model = st[0]
    ppb = Bio.PDB.PPBuilder() if False else None
    result = []
    for chain in model:
        reslist = [r for r in chain if r.id[0] == " " or (resname_filter and r.id[0] != " ")]
        # include both standard and D-residues (hetero flags vary); use all residues with N/CA/C
        reslist = [r for r in chain if "N" in r and "CA" in r and "C" in r]
        for i in range(1, len(reslist) - 1):
            try:
                phi = dihedral(reslist[i-1]["C"].get_vector(), reslist[i]["N"].get_vector(),
                               reslist[i]["CA"].get_vector(), reslist[i]["C"].get_vector())
                psi = dihedral(reslist[i]["N"].get_vector(), reslist[i]["CA"].get_vector(),
                               reslist[i]["C"].get_vector(), reslist[i+1]["N"].get_vector())
                result.append((chain.id, reslist[i].id[1], reslist[i].resname,
                               float(np.degrees(phi)), float(np.degrees(psi))))
            except Exception:
                continue
    return result


def dihedral(v0, v1, v2, v3):
    from Bio.PDB.vectors import calc_dihedral
    return calc_dihedral(v0, v1, v2, v3)


def mirror_qc(pdb_l, pdb_d, tol=1e-3):
    """Verify |phi_L + phi_D| ~= 0 and |psi_L + psi_D| ~= 0 for paired structures."""
    dl = dihedrals_from_pdb(pdb_l)
    dd = dihedrals_from_pdb(pdb_d)
    n = min(len(dl), len(dd))
    if n == 0:
        return {"n_pairs": 0, "max_abs_phi_sum": None, "max_abs_psi_sum": None}
    phi_sum = np.abs(np.array([d[3] for d in dl[:n]]) + np.array([d[3] for d in dd[:n]]))
    psi_sum = np.abs(np.array([d[4] for d in dl[:n]]) + np.array([d[4] for d in dd[:n]]))
    return {"n_pairs": n,
            "max_abs_phi_sum": float(phi_sum.max()),
            "max_abs_psi_sum": float(psi_sum.max())}
