import os
import csv
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from tqdm import tqdm

def calculate_pairwise_similarity(seq1, seq2):
    """
    Calculate the pairwise similarity between two sequences using PairwiseAligner.
    
    Parameters:
    seq1 (str): First DNA sequence.
    seq2 (str): Second DNA sequence.
    
    Returns:
    float: Similarity as a percentage.
    """
    aligner = PairwiseAligner()
    alignment = aligner.align(seq1, seq2)[0]  # Get the best alignment
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    length = max(len(seq1), len(seq2))
    similarity = matches / length * 100
    return similarity

def average_pairwise_similarity(sequences):
    """
    Calculate the average pairwise similarity for a list of sequences.
    
    Parameters:
    sequences (list): List of DNA sequences.
    
    Returns:
    float: Average pairwise similarity as a percentage.
    """
    if len(sequences) < 2:
        return 100.0 if sequences else 0.0

    total_similarity = 0.0
    count = 0

    with tqdm(total=len(sequences) * (len(sequences) - 1) // 2, desc="Calculating similarities", leave=False) as pbar:
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                similarity = calculate_pairwise_similarity(sequences[i], sequences[j])
                total_similarity += similarity
                count += 1
                pbar.update(1)

    average_similarity = total_similarity / count if count > 0 else 0.0
    return average_similarity

def process_fasta_files_in_directory(directory):
    """
    Process all FASTA files in a directory to calculate average pairwise similarity.
    
    Parameters:
    directory (str): Path to the directory containing FASTA files.
    
    Returns:
    list: A list of tuples containing filename and average pairwise similarity.
    """
    results = []
    files = [f for f in os.listdir(directory) if f.endswith(".fasta") or f.endswith(".fa")]
    for filename in tqdm(files, desc="Processing files"):
        filepath = os.path.join(directory, filename)
        sequences = [str(record.seq) for record in SeqIO.parse(filepath, "fasta")]
        avg_similarity = average_pairwise_similarity(sequences)
        results.append((filename, avg_similarity))
    return results

def save_results_to_csv(results, output_file):
    """
    Save the average pairwise similarity results to a CSV file.
    
    Parameters:
    results (list): A list of tuples containing filename and average pairwise similarity.
    output_file (str): Path to the output CSV file.
    """
    with open(output_file, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Filename", "Average Pairwise Similarity (%)"])
        for row in results:
            writer.writerow(row)

# Directory containing FASTA files
directory = '.'

# Output CSV file
output_file = 'average_pairwise_similarity_results.csv'

# Process FASTA files and save results to CSV
results = process_fasta_files_in_directory(directory)
save_results_to_csv(results, output_file)

print(f"Average pairwise similarity results have been saved to {output_file}")
