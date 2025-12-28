# Stub RL placer: reads devices CSV and emits x positions spaced evenly.
# Usage: python -m rl_placer.main --input in.csv --output out.csv

import argparse
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    devices = []
    with open(args.input, newline='') as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            devices.append(row)

    # Simple policy: x = index * pitch + margin
    pitch = 0.19  # um
    margin = 0.19
    with open(args.output, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['id', 'fold', 'x'])
        for i, d in enumerate(devices):
            w.writerow([d['id'], d['fold'], margin + i * pitch])


if __name__ == '__main__':
    main()
