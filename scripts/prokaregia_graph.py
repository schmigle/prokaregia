import os, argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

# Set directories
input_dir = 'checkm2_combined'
output_dir = 'checkm2_graphs'
os.makedirs(output_dir, exist_ok=True)

# Define color map and palette
color_map = {'polished': 'Polished', 'dastool': 'Unpolished', 'rosella': 'Rosella', 'metabat':'Metabat', 'maxbin':'Maxbin', 'semibin': 'Semibin', 'metacoag': 'Metacoag'}
color_palette = {'Polished': '#ABABAB', 'Unpolished': '#79AAE6', 'Rosella': '#5C8EAC', 'Metabat': '#D4CFCF', 'Maxbin': '#A3D977', 'Semibin': '#F2E6B6', 'Metacoag': '#E4A8A1'}

# Extract group from filename
def extract_group(filename):
    parts = filename.split(".tsv")[0]

# Group files by prefix
files_grouped = defaultdict(list)
for file in os.listdir(input_dir):
    if file.endswith('.tsv'):
        files_grouped[extract_group(file)].append(file)

# Plotting function
def plot_data(data_list, y_column, title, output_filename, y_min, y_max, y_ticks):
    plt.figure(figsize=(8.5, 8))

    # Determine if any line has 5 or more points
    any_line_with_5_or_more = any(len(data['Data']) >= 5 for data in data_list)

    for data in data_list:
        label_key = next((key for key in color_map if key in data['File']), data['File'])
        label = color_map.get(label_key, label_key)
        color = color_palette.get(label)
        linestyle = '-' if any_line_with_5_or_more else ('o-' if len(data['Data']) < 5 else '-')
        plt.plot(data['Data'].index + 1, data['Data'][y_column], linestyle, \
                 color=color, label=label, linewidth=1.5)

    plt.title(title, fontsize=24)
    plt.xlabel("Rank", fontsize=18)
    plt.ylabel(f"{y_column} (%)", fontsize=18)
    plt.ylim(y_min, y_max)
    plt.yticks(y_ticks, fontsize=12)

    max_rank = max(len(data['Data']) for data in data_list)
    gap = max(1, (max_rank - 1) // 10 + 1)
    x_ticks = range(gap, max_rank + 1, gap)
    plt.xticks(x_ticks, fontsize=14)
    plt.legend(bbox_to_anchor=(0.5, -0.27), loc='lower center', fontsize=16, ncol=4)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_filename))
    plt.close()

# Process and plot for each group
for group, files in files_grouped.items():
    completeness_data, contamination_data = [], []
    for file in files:
        df = pd.read_csv(os.path.join(input_dir, file), sep='\t')
        df = df[(df['Contamination'] <= 10) & (df['Completeness'] >= 50)]
        if not df.empty:
            completeness_data.append({'File': file, 'Data': df.sort_values(by='Completeness', ascending=False).reset_index()})
            contamination_data.append({'File': file, 'Data': df.sort_values(by='Contamination').reset_index()})
    if completeness_data:
        plot_data(completeness_data, 'Completeness', f"Completeness", f"Completeness.png", 48, 102, range(50, 101, 10))
    if contamination_data:
        plot_data(contamination_data, 'Contamination', f"Contamination", f"Contamination.png", -0.5, 11, range(0, 11, 2))
