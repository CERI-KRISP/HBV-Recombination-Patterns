import os
import csv
from Bio import SeqIO

def calculate_similarity_per_position(sequences):
    """
    Calculate the sequence similarity for every position across the genome.
    
    Parameters:
    sequences (list): List of DNA sequences.
    
    Returns:
    dict: A dictionary with position as key and similarity as value.
    """
    position_similarity = {}
    max_length = max(len(seq) for seq in sequences)
    
    for pos in range(max_length):
        base_counts = {}
        valid_bases = 0
        
        for seq in sequences:
            if pos < len(seq):
                base = seq[pos]
                if base in ['G', 'C', 'A', 'T']:
                    valid_bases += 1
                    if base in base_counts:
                        base_counts[base] += 1
                    else:
                        base_counts[base] = 1
        
        if valid_bases > 0:
            most_common_base_count = max(base_counts.values())
            similarity = most_common_base_count / valid_bases * 100
        else:
            similarity = 0.0
        
        position_similarity[f'Position_{pos + 1}'] = similarity
    
    return position_similarity

def process_fasta_files_in_directory(directory):
    """
    Process all FASTA files in a directory to calculate sequence similarity per position.
    
    Parameters:
    directory (str): Path to the directory containing FASTA files.
    
    Returns:
    list: A list of dictionaries containing filename and sequence similarity per position.
    """
    results = []
    for filename in os.listdir(directory):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            filepath = os.path.join(directory, filename)
            sequences = [str(record.seq) for record in SeqIO.parse(filepath, "fasta")]
            similarity_per_position = calculate_similarity_per_position(sequences)
            result = {'Filename': filename, **similarity_per_position}
            results.append(result)
    return results

def save_results_to_csv(results, output_file):
    """
    Save the sequence similarity results to a CSV file.
    
    Parameters:
    results (list): A list of dictionaries containing filename and sequence similarity per position.
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
output_file = 'sequence_similarity_per_position_results.csv'

# Process FASTA files and save results to CSV
results = process_fasta_files_in_directory(directory)
save_results_to_csv(results, output_file)

print(f"Sequence similarity per position results have been saved to {output_file}")