"""CPU interface descriptors for protein-peptide complexes.

Lightweight, dependency-light descriptors computed from the co-crystal structure:
contact counts, interface area estimate, H-bonds, salt bridges, charge
complementarity, hydrophobic fraction. These feed the Mirror-Score calibration
alongside structure-prediction confidence metrics.
"""
import numpy as np
from Bio.PDB import PDBParser
import warnings
warnings.filterwarnings("ignore")

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}
POSITIVE = {"ARG", "LYS", "HIS"}
NEGATIVE = {"ASP", "GLU"}
POLAR = {"SER", "THR", "ASN", "GLN", "CYS", "TYR", "HIS"}

# D-residue CCD codes -> L-equivalent for composition analysis
D_TO_L_MAP = {
    "DAL": "ALA", "DAR": "ARG", "DSG": "ASN", "DAS": "ASP", "DCY": "CYS",
    "DGN": "GLN", "DGL": "GLU", "DHI": "HIS", "DIL": "ILE", "DLE": "LEU",
    "DLY": "LYS", "MED": "MET", "DPN": "PHE", "DPR": "PRO", "DSN": "SER",
    "DTH": "THR", "DTR": "TRP", "DTY": "TYR", "DVA": "VAL", "DSP": "ASP",
}


def _resname(res):
    rn = res.resname.strip().upper()
    return D_TO_L_MAP.get(rn, rn)


def _cb_coord(res):
    if "CB" in res:
        return res["CB"].coord
    if _resname(res) == "GLY" and "CA" in res:
        return res["CA"].coord
    return res["CA"].coord if "CA" in res else None


def load_chain_residues(pdb_path, chain_id=None):
    parser = PDBParser(QUIET=True)
    st = parser.get_structure("x", pdb_path)
    model = st[0]
    chains = {}
    for ch in model:
        reslist = [r for r in ch if "CA" in r]
        if reslist:
            chains[ch.id] = reslist
    return chains


def interface_descriptors(pdb_path, peptide_chain, target_chains=None, cutoff=5.0):
    """Compute interface descriptors between a peptide chain and target chain(s).

    Works for both L-peptide/L-protein and D-peptide/L-protein complexes
    (D-residues appear as HETATM/alt CCD codes but have N/CA/C/CB atoms).
    """
    chains = load_chain_residues(pdb_path)
    if peptide_chain not in chains:
        raise KeyError(f"chain {peptide_chain} not found; have {list(chains)}")
    pep = chains[peptide_chain]
    if target_chains is None:
        target_chains = [c for c in chains if c != peptide_chain]
    tgt = [r for c in target_chains if c in chains for r in chains[c]]

    pep_atoms = [a for r in pep for a in r if a.element != "H"]
    tgt_atoms = [a for r in tgt for a in r if a.element != "H"]
    if not pep_atoms or not tgt_atoms:
        return None

    pep_xyz = np.array([a.coord for a in pep_atoms]).reshape(-1, 3)
    tgt_xyz = np.array([a.coord for a in tgt_atoms]).reshape(-1, 3)
    d = np.linalg.norm(pep_xyz[:, None, :] - tgt_xyz[None, :, :], axis=2)
    n_contacts = int((d < cutoff).sum())

    # per-residue contacts
    pep_res_contacts = set()
    tgt_res_contacts = set()
    idx = 0
    pep_atom_res = []
    for r in pep:
        for a in r:
            if a.element != "H":
                pep_atom_res.append(id(r))
                idx += 1
    contact_mask = d < cutoff
    pep_hit = contact_mask.any(axis=1)
    tgt_hit = contact_mask.any(axis=0)
    # map atom indices back to residues
    pep_res_ids = []
    for r in pep:
        n_at = sum(1 for a in r if a.element != "H")
        pep_res_ids.extend([r.id[1]] * n_at)
    tgt_res_ids = []
    for r in tgt:
        n_at = sum(1 for a in r if a.element != "H")
        tgt_res_ids.extend([r.id[1]] * n_at)
    pep_res_contacts = len(set(np.array(pep_res_ids)[pep_hit]))
    tgt_res_contacts = len(set(np.array(tgt_res_ids)[tgt_hit]))

    # H-bonds (donor-acceptor heavy-atom distance < 3.5 A, simplified geometric criterion)
    hb = 0
    for i, a in enumerate(pep_atoms):
        if a.element in ("N", "O"):
            close = np.where(d[i] < 3.5)[0]
            for j in close:
                b = tgt_atoms[j]
                if b.element in ("N", "O") and (a.name.startswith("N") and b.name.startswith("O")
                                                or a.name.startswith("O") and b.name.startswith("N")):
                    hb += 1

    # composition-based descriptors
    pep_types = [_resname(r) for r in pep]
    tgt_types = [_resname(r) for r in tgt]
    pep_charge = sum(1 for t in pep_types if t in POSITIVE) - sum(1 for t in pep_types if t in NEGATIVE)
    tgt_charge = sum(1 for t in tgt_types if t in POSITIVE) - sum(1 for t in tgt_types if t in NEGATIVE)
    pep_hydro = sum(1 for t in pep_types if t in HYDROPHOBIC) / max(len(pep_types), 1)

    # salt bridges: opposite-charge residue pairs with any atom pair < 4.5 A
    salt = 0
    pos_res = [r for r in pep if _resname(r) in POSITIVE]
    neg_res = [r for r in tgt if _resname(r) in NEGATIVE]
    for rp in pos_res:
        for rn in neg_res:
            pa = np.array([a.coord for a in rp if a.element != "H"]).reshape(-1, 3)
            na = np.array([a.coord for a in rn if a.element != "H"]).reshape(-1, 3)
            if pa.size and na.size and np.linalg.norm(pa[:, None] - na[None, :], axis=2).min() < 4.5:
                salt += 1
    # and the reverse (peptide negative, target positive)
    neg_res_p = [r for r in pep if _resname(r) in NEGATIVE]
    pos_res_t = [r for r in tgt if _resname(r) in POSITIVE]
    for rp in neg_res_p:
        for rn in pos_res_t:
            pa = np.array([a.coord for a in rp if a.element != "H"]).reshape(-1, 3)
            na = np.array([a.coord for a in rn if a.element != "H"]).reshape(-1, 3)
            if pa.size and na.size and np.linalg.norm(pa[:, None] - na[None, :], axis=2).min() < 4.5:
                salt += 1

    return {
        "n_atom_contacts": n_contacts,
        "n_pep_res_in_contact": pep_res_contacts,
        "n_tgt_res_in_contact": tgt_res_contacts,
        "n_hbonds": hb,
        "n_salt_bridges": salt,
        "peptide_net_charge": pep_charge,
        "target_net_charge": tgt_charge,
        "charge_complementarity": tgt_charge - pep_charge,
        "peptide_hydrophobic_frac": pep_hydro,
        "peptide_len": len(pep),
    }
