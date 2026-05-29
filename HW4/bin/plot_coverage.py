#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

depth_file = sys.argv[1]
output_file = sys.argv[2]

positions = []
depths = []

with open(depth_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        positions.append(int(parts[1]))
        depths.append(int(parts[2]))

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(positions, depths, linewidth=0.5, color='steelblue')
ax.fill_between(positions, depths, alpha=0.3, color='steelblue')
ax.set_xlabel('Position')
ax.set_ylabel('Coverage depth')
ax.set_title('Coverage')

plt.tight_layout()
plt.savefig(output_file, dpi=150)
