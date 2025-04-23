import os
import sys
from collections import defaultdict

def summarize_density(directory):
    contig_counts = defaultdict(int)
    all_contigs = set()

    if not os.path.isdir(directory):
        print(f"Error: '{directory}' is not a valid directory", file=sys.stderr)
        sys.exit(1)

    for filename in os.listdir(directory):
        if filename.endswith('.stats'):
            filepath = os.path.join(directory, filename)
            
            try:
                with open(filepath, 'r') as file:
                    first_line = file.readline().strip()
                    fields = first_line.split('\t')
                    
                    if len(fields) >= 5:
                        line = file.readline().strip()
                        while line:
                            fields = line.split('\t')
                            contig_name = fields[0]
                            all_contigs.add(contig_name)
                            density = float(fields[4])
                            if density > 0.1:
                                contig_counts[contig_name] += 1
                            line = file.readline().strip()
            except (IOError, ValueError, IndexError) as e:
                print(f"Skipping {filename} due to error: {e}", file=sys.stderr)

    print("Contig Name\tFiles with Density > 0.1")
    print("------------------------------------")
    max_count = 0
    for contig in sorted(all_contigs):
        if contig_counts.get(contig, 0) > max_count:
            max_count = contig_counts.get(contig, 0)
        print(f"{contig}\t{contig_counts.get(contig, 0)}")

    print("Max count:", max_count)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <directory_path>", file=sys.stderr)
        sys.exit(1)
    
    summarize_density(sys.argv[1])