"""Process completed Boltz-2 batch results: extract mirror-space features
for every prediction and merge into data/benchmark/boltz_features.json.

Usage: python process_boltz_batch.py <batch_outputs_root>
Expects <root>/outputs/boltz_results_<NAME>/predictions/<NAME>/...
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mirrorscore.boltz_features import extract_features  # noqa: E402


def main(root):
    out_path = os.path.join(os.path.dirname(root), 'boltz_features.json')
    # prefer the canonical location
    out_path = '/workspace/mirrorscore/data/benchmark/boltz_features.json'
    feats_all = json.load(open(out_path)) if os.path.exists(out_path) else {}

    outputs = os.path.join(root, 'outputs')
    runs = [d for d in os.listdir(outputs) if d.startswith('boltz_results_')]
    print(f'found {len(runs)} boltz result dirs')

    ok, fail = 0, []
    for d in sorted(runs):
        name = d.replace('boltz_results_', '')
        pred = os.path.join(outputs, d, 'predictions', name)
        if not os.path.isdir(pred):
            fail.append((name, 'no predictions dir'))
            continue
        try:
            feats_all[name] = extract_features(pred, name)
            ok += 1
        except Exception as e:
            fail.append((name, str(e)[:120]))

    json.dump(feats_all, open(out_path, 'w'), indent=1)
    print(f'extracted {ok} ok, {len(fail)} failed')
    for name, err in fail:
        print(' FAIL', name, err)


if __name__ == '__main__':
    main(sys.argv[1])
