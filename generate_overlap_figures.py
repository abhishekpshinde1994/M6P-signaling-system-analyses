#!/usr/bin/env python3

"""
Generate overlap analysis figures for co-IP proteomics data
Creates three Venn diagrams with corresponding heatmaps:
A. Overall protein overlap between GNPT, MPR3, and MPR7
B. Lysosomal protein overlap
C. Secretome protein overlap
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Circle, Ellipse
import numpy as np
import seaborn as sns

# Load data
df = pd.read_csv('/mnt/user-data/outputs/interaction_summary_with_annotations.csv')

# Extract protein sets for each bait
gnpt_proteins = set(df[df['baits_interacting'].str.contains('GNPT', na=False)]['protein_id'])
mpr3_proteins = set(df[df['baits_interacting'].str.contains('MPR3', na=False)]['protein_id'])
mpr7_proteins = set(df[df['baits_interacting'].str.contains('MPR7', na=False)]['protein_id'])

# Calculate all overlaps for main figure
all_three = gnpt_proteins & mpr3_proteins & mpr7_proteins
gnpt_mpr3_only = (gnpt_proteins & mpr3_proteins) - mpr7_proteins
gnpt_mpr7_only = (gnpt_proteins & mpr7_proteins) - mpr3_proteins
mpr3_mpr7_only = (mpr3_proteins & mpr7_proteins) - gnpt_proteins
gnpt_only = gnpt_proteins - mpr3_proteins - mpr7_proteins
mpr3_only = mpr3_proteins - gnpt_proteins - mpr7_proteins
mpr7_only = mpr7_proteins - gnpt_proteins - mpr3_proteins

# Lysosome overlaps
gnpt_lyso = set(df[(df['baits_interacting'].str.contains('GNPT', na=False)) & 
                   (df['lysosome'] == '+')]['protein_id'])
mpr3_lyso = set(df[(df['baits_interacting'].str.contains('MPR3', na=False)) & 
                   (df['lysosome'] == '+')]['protein_id'])
mpr7_lyso = set(df[(df['baits_interacting'].str.contains('MPR7', na=False)) & 
                   (df['lysosome'] == '+')]['protein_id'])

all_three_lyso = gnpt_lyso & mpr3_lyso & mpr7_lyso
gnpt_mpr3_only_lyso = (gnpt_lyso & mpr3_lyso) - mpr7_lyso
gnpt_mpr7_only_lyso = (gnpt_lyso & mpr7_lyso) - mpr3_lyso
mpr3_mpr7_only_lyso = (mpr3_lyso & mpr7_lyso) - gnpt_lyso
gnpt_only_lyso = gnpt_lyso - mpr3_lyso - mpr7_lyso
mpr3_only_lyso = mpr3_lyso - gnpt_lyso - mpr7_lyso
mpr7_only_lyso = mpr7_lyso - gnpt_lyso - mpr3_lyso

# Secretome overlaps
gnpt_sec = set(df[(df['baits_interacting'].str.contains('GNPT', na=False)) & 
                  (df['secretome'] == '+')]['protein_id'])
mpr3_sec = set(df[(df['baits_interacting'].str.contains('MPR3', na=False)) & 
                  (df['secretome'] == '+')]['protein_id'])
mpr7_sec = set(df[(df['baits_interacting'].str.contains('MPR7', na=False)) & 
                  (df['secretome'] == '+')]['protein_id'])

all_three_sec = gnpt_sec & mpr3_sec & mpr7_sec
gnpt_mpr3_only_sec = (gnpt_sec & mpr3_sec) - mpr7_sec
gnpt_mpr7_only_sec = (gnpt_sec & mpr7_sec) - mpr3_sec
mpr3_mpr7_only_sec = (mpr3_sec & mpr7_sec) - gnpt_sec
gnpt_only_sec = gnpt_sec - mpr3_sec - mpr7_sec
mpr3_only_sec = mpr3_sec - gnpt_sec - mpr7_sec
mpr7_only_sec = mpr7_sec - gnpt_sec - mpr3_sec


def plot_venn3_with_heatmap(counts_dict, title, ax_venn, ax_heat, dataset_name=""):
    """
    Plot a 3-way Venn diagram with corresponding heatmap
    """
    # Unpack counts
    gnpt_only = counts_dict['gnpt_only']
    mpr3_only = counts_dict['mpr3_only']
    mpr7_only = counts_dict['mpr7_only']
    gnpt_mpr3 = counts_dict['gnpt_mpr3']
    gnpt_mpr7 = counts_dict['gnpt_mpr7']
    mpr3_mpr7 = counts_dict['mpr3_mpr7']
    all_three = counts_dict['all_three']
    
    # Draw Venn diagram
    # Circle parameters (similar to your reference)
    circle_gnpt = Circle((0.35, 0.5), 0.25, alpha=0.5, color='lightblue', ec='black', linewidth=1.5)
    circle_mpr3 = Circle((0.55, 0.6), 0.25, alpha=0.5, color='lightgreen', ec='black', linewidth=1.5)
    circle_mpr7 = Circle((0.55, 0.4), 0.25, alpha=0.5, color='lightcoral', ec='black', linewidth=1.5)
    
    ax_venn.add_patch(circle_gnpt)
    ax_venn.add_patch(circle_mpr3)
    ax_venn.add_patch(circle_mpr7)
    
    # Add numbers
    fontsize = 14
    fontweight = 'bold'
    
    # Only regions
    ax_venn.text(0.20, 0.5, str(gnpt_only), ha='center', va='center', 
                fontsize=fontsize, fontweight=fontweight)
    ax_venn.text(0.60, 0.73, str(mpr3_only), ha='center', va='center', 
                fontsize=fontsize, fontweight=fontweight)
    ax_venn.text(0.60, 0.27, str(mpr7_only), ha='center', va='center', 
                fontsize=fontsize, fontweight=fontweight)
    
    # Two-way overlaps
    ax_venn.text(0.43, 0.58, str(gnpt_mpr3), ha='center', va='center', 
                fontsize=fontsize, fontweight=fontweight)
    ax_venn.text(0.43, 0.42, str(gnpt_mpr7), ha='center', va='center', 
                fontsize=fontsize, fontweight=fontweight)
    ax_venn.text(0.65, 0.5, str(mpr3_mpr7), ha='center', va='center', 
                fontsize=fontsize, fontweight=fontweight)
    
    # Three-way overlap
    ax_venn.text(0.50, 0.5, str(all_three), ha='center', va='center', 
                fontsize=fontsize, fontweight=fontweight)
    
    # Labels
    ax_venn.text(0.15, 0.75, 'GNPT', ha='center', va='center', 
                fontsize=12, fontweight='bold')
    ax_venn.text(0.75, 0.75, 'MPR3', ha='center', va='center', 
                fontsize=12, fontweight='bold')
    ax_venn.text(0.75, 0.25, 'MPR7', ha='center', va='center', 
                fontsize=12, fontweight='bold')
    
    ax_venn.set_xlim(0, 1)
    ax_venn.set_ylim(0, 1)
    ax_venn.set_aspect('equal')
    ax_venn.axis('off')
    
    # Create heatmap data matching the reference figure format
    # The heatmap shows: how many proteins from each row dataset overlap with column dataset
    # Rows: MPR7, MPR3, GNPT (top to bottom)
    # Columns: GNPT, MPR3, MPR7 (left to right)
    
    # Calculate totals for each bait
    gnpt_total = gnpt_only + gnpt_mpr3 + gnpt_mpr7 + all_three
    mpr3_total = mpr3_only + gnpt_mpr3 + mpr3_mpr7 + all_three
    mpr7_total = mpr7_only + gnpt_mpr7 + mpr3_mpr7 + all_three
    
    # Build matrix: rows = MPR7, MPR3, GNPT; cols = GNPT, MPR3, MPR7
    heatmap_data = np.array([
        [gnpt_mpr7 + all_three,  # MPR7 proteins that overlap with GNPT
         mpr3_mpr7 + all_three,  # MPR7 proteins that overlap with MPR3
         mpr7_total],             # All MPR7 proteins
        
        [gnpt_mpr3 + all_three,  # MPR3 proteins that overlap with GNPT
         mpr3_total,             # All MPR3 proteins
         mpr3_mpr7 + all_three], # MPR3 proteins that overlap with MPR7
        
        [gnpt_total,             # All GNPT proteins
         gnpt_mpr3 + all_three,  # GNPT proteins that overlap with MPR3
         gnpt_mpr7 + all_three]  # GNPT proteins that overlap with MPR7
    ])
    
    # Create heatmap
    sns.heatmap(heatmap_data, annot=True, fmt='d', cmap='RdYlBu_r', 
                cbar_kws={'label': 'Protein Count'},
                xticklabels=['GNPT', 'MPR3', 'MPR7'],
                yticklabels=['MPR7', 'MPR3', 'GNPT'],
                ax=ax_heat, linewidths=1.5, linecolor='white',
                square=True, cbar=True,
                annot_kws={'fontsize': 12, 'fontweight': 'bold'})
    
    ax_heat.set_xlabel('Dataset', fontsize=11, fontweight='bold')
    ax_heat.set_ylabel('Dataset', fontsize=11, fontweight='bold')
    ax_heat.set_xticklabels(ax_heat.get_xticklabels(), rotation=0, fontsize=10)
    ax_heat.set_yticklabels(ax_heat.get_yticklabels(), rotation=0, fontsize=10)
    
    # Add title above heatmap
    ax_heat.set_title(dataset_name, fontsize=10, fontweight='bold', pad=12)
    
    return ax_venn, ax_heat


# Create the figure
fig = plt.figure(figsize=(16, 14))

# Panel A - Overall overlap
ax_a_venn = plt.subplot(3, 2, 1)
ax_a_heat = plt.subplot(3, 2, 2)

counts_a = {
    'gnpt_only': len(gnpt_only),
    'mpr3_only': len(mpr3_only),
    'mpr7_only': len(mpr7_only),
    'gnpt_mpr3': len(gnpt_mpr3_only),
    'gnpt_mpr7': len(gnpt_mpr7_only),
    'mpr3_mpr7': len(mpr3_mpr7_only),
    'all_three': len(all_three)
}

plot_venn3_with_heatmap(counts_a, 'A. CoIP Overlap analysis', 
                        ax_a_venn, ax_a_heat, 
                        dataset_name='coIP Protein Dataset Relationships')

# Panel B - Lysosome overlap
ax_b_venn = plt.subplot(3, 2, 3)
ax_b_heat = plt.subplot(3, 2, 4)

counts_b = {
    'gnpt_only': len(gnpt_only_lyso),
    'mpr3_only': len(mpr3_only_lyso),
    'mpr7_only': len(mpr7_only_lyso),
    'gnpt_mpr3': len(gnpt_mpr3_only_lyso),
    'gnpt_mpr7': len(gnpt_mpr7_only_lyso),
    'mpr3_mpr7': len(mpr3_mpr7_only_lyso),
    'all_three': len(all_three_lyso)
}

plot_venn3_with_heatmap(counts_b, 'B. Lysosomal proteins in coIP Overlap analysis', 
                        ax_b_venn, ax_b_heat,
                        dataset_name='Lysosome-Overlapping Protein Relationships')

# Panel C - Secretome overlap
ax_c_venn = plt.subplot(3, 2, 5)
ax_c_heat = plt.subplot(3, 2, 6)

counts_c = {
    'gnpt_only': len(gnpt_only_sec),
    'mpr3_only': len(mpr3_only_sec),
    'mpr7_only': len(mpr7_only_sec),
    'gnpt_mpr3': len(gnpt_mpr3_only_sec),
    'gnpt_mpr7': len(gnpt_mpr7_only_sec),
    'mpr3_mpr7': len(mpr3_mpr7_only_sec),
    'all_three': len(all_three_sec)
}

plot_venn3_with_heatmap(counts_c, 'C. Secreted proteins in coIP Overlap analysis', 
                        ax_c_venn, ax_c_heat,
                        dataset_name='Secretome-Overlapping Protein Relationships')

# Add panel labels
fig.text(0.02, 0.95, 'A. CoIP Overlap analysis', fontsize=14, fontweight='bold')
fig.text(0.02, 0.63, 'B. Lysosomal proteins in coIP Overlap analysis', 
         fontsize=14, fontweight='bold')
fig.text(0.02, 0.31, 'C. Secreted proteins in coIP Overlap analysis', 
         fontsize=14, fontweight='bold')

plt.tight_layout(rect=[0, 0, 1, 0.98])
plt.savefig('/mnt/user-data/outputs/FigureS4_overlap_analysis.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/mnt/user-data/outputs/FigureS4_overlap_analysis.png', dpi=300, bbox_inches='tight')

print("Figure saved successfully!")
print("\n=== Summary ===")
print("A. Overall CoIP overlaps:")
print(f"  GNPT only: {len(gnpt_only)}")
print(f"  MPR3 only: {len(mpr3_only)}")
print(f"  MPR7 only: {len(mpr7_only)}")
print(f"  GNPT-MPR3: {len(gnpt_mpr3_only)}")
print(f"  GNPT-MPR7: {len(gnpt_mpr7_only)}")
print(f"  MPR3-MPR7: {len(mpr3_mpr7_only)}")
print(f"  All three: {len(all_three)}")

print("\nB. Lysosome overlaps:")
print(f"  Total: {len(gnpt_only_lyso) + len(mpr3_only_lyso) + len(mpr7_only_lyso) + len(gnpt_mpr3_only_lyso) + len(gnpt_mpr7_only_lyso) + len(mpr3_mpr7_only_lyso) + len(all_three_lyso)}")

print("\nC. Secretome overlaps:")
print(f"  Total: {len(gnpt_only_sec) + len(mpr3_only_sec) + len(mpr7_only_sec) + len(gnpt_mpr3_only_sec) + len(gnpt_mpr7_only_sec) + len(mpr3_mpr7_only_sec) + len(all_three_sec)}")
