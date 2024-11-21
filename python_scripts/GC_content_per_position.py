import os
import csv
from Bio import SeqIO

def calculate_gc_content_per_position(sequences):
    """
    Calculate the GC content for every position across the genome.
    
    Parameters:
    sequences (list): List of DNA sequences.
    
    Returns:
    dict: A dictionary with position as key and GC content as value.
    """
    position_gc_content = {}
    max_length = max(len(seq) for seq in sequences)
    
    for pos in range(max_length):
        g_count = 0
        c_count = 0
        valid_bases = 0
        
        for seq in sequences:
            if pos < len(seq):
                base = seq[pos]
                if base == 'G':
                    g_count += 1
                elif base == 'C':
                    c_count += 1
                if base in ['G', 'C', 'A', 'T']:
                    valid_bases += 1
        
        if valid_bases > 0:
            gc_content = (g_count + c_count) / valid_bases * 100
        else:
            gc_content = 0.0
        
        position_gc_content[f'Position_{pos + 1}'] = gc_content
    
    return position_gc_content

def process_fasta_files_in_directory(directory):
    """
    Process all FASTA files in a directory to calculate GC content per position.
    
    Parameters:
    directory (str): Path to the directory containing FASTA files.
    
    Returns:
    list: A list of dictionaries containing filename and GC content per position.
    """
    results = []
    for filename in os.listdir(directory):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            filepath = os.path.join(directory, filename)
            sequences = [str(record.seq) for record in SeqIO.parse(filepath, "fasta")]
            gc_content_per_position = calculate_gc_content_per_position(sequences)
            result = {'Filename': filename, **gc_content_per_position}
            results.append(result)
    return results

def save_results_to_csv(results, output_file):
    """
    Save the GC content results to a CSV file.
    
    Parameters:
    results (list): A list of dictionaries containing filename and GC content per position.
    output_file (str): Path to the output CSV file.
    """
    # Determine the maximum number of positions to set the header
    max_positions = max(len(result) - 1 for result in results)  # -1 to exclude the Filename key
    fieldnames = ['Filename'] + [f'Position_{i + 1}' for i in range(max_positions)]
    
    with open(output_file, mode='w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            # Fill in missing positions with empty strings for consistency
            for i in range(1, max_positions + 1):
                if f'Position_{i}' not in result:
                    result[f'Position_{i}'] = ''
            writer.writerow(result)

# Directory containing FASTA files
directory = '.'

# Output CSV file
output_file = 'gc_content_per_position_results.csv'

# Process FASTA files and save results to CSV
results = process_fasta_files_in_directory(directory)
save_results_to_csv(results, output_file)

print(f"GC content per position results have been saved to {output_file}")
