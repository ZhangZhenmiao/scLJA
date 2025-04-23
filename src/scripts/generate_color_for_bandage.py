import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def normalize_and_map_colors(df):
    # Normalize the density values using min-max scaling
    # nonzero = density_values[density_values > 0]
    # vmin = nonzero.min() if not nonzero.empty else 0
    # vmax = nonzero.max() if not nonzero.empty else 1
    
    # norm = Normalize(vmin=vmin, vmax=vmax)
    # cmap = plt.cm.Reds  # Use the "Reds" colormap

    colors = []
    for i in range(len(df)):
        if df["density"][i] < 0.1:
            colors.append("#cccacb")  # gray for zero
        # else:
        #     norm_val = (val - vmin) / (vmax - vmin) if vmax > vmin else 1.0
        #     color = ScalarMappable(norm=norm, cmap=cmap).to_rgba(val, bytes=True)
        #     colors.append('#{:02x}{:02x}{:02x}'.format(color[0], color[1], color[2]))  # RGB hex
        else:
            colors.append("#ff0000")
    return colors


def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_tsv> <output_tsv>")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_txt = sys.argv[2]

    df = pd.read_csv(input_csv, sep='\t', dtype=str)

    if 'contig_name' not in df.columns or 'density' not in df.columns:
        print("Input file must contain 'contig_name' and 'density' columns.")
        sys.exit(1)

    df['density'] = df['density'].astype(float)
    df['read_count'] = df['read_count'].astype(int)

    df_output = pd.DataFrame()
    df_output['contig_name'] = df['contig_name']  # keep it as string without type conversion
    df_output['color'] = normalize_and_map_colors(df)

    # Write to tab-separated file with header
    df_output.to_csv(output_txt, index=False, header=True, sep='\t')


if __name__ == "__main__":
    main()
