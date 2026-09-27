"""Extract mirror-space Boltz-2 confidence features for a complex prediction.

Features (all computed in mirror space, where both chains are L-form):
- iptm: interface TM-score (global interface confidence)
- pair_chains_iptm_min: min of cross-chain ipTM values (direction-averaged)
- complex_iplddt: interface pLDDT
- complex_pde / complex_ipde: predicted distance error (global / interface)
- if_pae: mean PAE between peptide and target tokens (interface PAE)
- if_pae_pep: mean PAE of peptide-internal tokens (peptide self-structure)
- confidence_score: Boltz aggregate
"""
import json
import numpy as np


def _chain_token_ranges(cif_path):
    """Parse token (residue) ranges per chain from the predicted CIF.

    Uses label_asym_id ordering, which matches Boltz token indexing for
    polymer chains. Returns {chain_label: (start_token, n_tokens)}.
    """
    import gemmi
    doc = gemmi.read_structure(cif_path)
    model = doc[0]
    ranges = {}
    start = 0
    for chain in model:
        n = sum(1 for res in chain if not res.name.startswith('HOH'))
        ranges[chain.name] = (start, n)
        start += n
    return ranges


def extract_features(pred_dir, name):
    """pred_dir contains {name}_model_0.cif, confidence_{name}_model_0.json,
    pae_{name}_model_0.npz."""
    conf = json.load(open(f'{pred_dir}/confidence_{name}_model_0.json'))
    pae = np.load(f'{pred_dir}/pae_{name}_model_0.npz')['pae']

    feats = {
        'boltz_confidence': conf['confidence_score'],
        'boltz_ptm': conf['ptm'],
        'boltz_iptm': conf['iptm'],
        'boltz_complex_iplddt': conf['complex_iplddt'],
        'boltz_complex_pde': conf['complex_pde'],
        'boltz_complex_ipde': conf['complex_ipde'],
        'boltz_pair_chains_iptm_min': min(conf['pair_chains_iptm']['0'].values()),
    }

    # chain token ranges from CIF (chain 0 = target, chain 1 = peptide by YAML order)
    try:
        ranges = _chain_token_ranges(f'{pred_dir}/{name}_model_0.cif')
        keys = sorted(ranges, key=lambda k: ranges[k][0])
        tgt = ranges[keys[0]]
        pep = ranges[keys[1]]
        n_tgt, n_pep = tgt[1], pep[1]
        pep_slice = slice(tgt[1], tgt[1] + n_pep)
        tgt_slice = slice(0, n_tgt)

        if_pae = pae[pep_slice, tgt_slice].mean()
        pep_self_pae = pae[pep_slice, pep_slice].mean()
        feats['boltz_if_pae'] = float(if_pae)
        feats['boltz_pep_self_pae'] = float(pep_self_pae)
        feats['n_tgt_tokens'] = n_tgt
        feats['n_pep_tokens'] = n_pep
    except Exception as e:  # CIF parsing optional
        feats['boltz_if_pae'] = float('nan')
        feats['extract_error'] = str(e)

    return feats


if __name__ == '__main__':
    import sys
    feats = extract_features(sys.argv[1], sys.argv[2])
    print(json.dumps(feats, indent=1))
