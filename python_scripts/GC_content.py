import os
import csv
from Bio import SeqIO

def calculate_gc_content(sequence):
    """
    Calculate the GC content of a DNA sequence.
    
    Parameters:
    sequence (str): DNA sequence.
    
    Returns:
    float: GC content as a percentage.
    """
    g = sequence.count('G')
    c = sequence.count('C')
    gc_content = (g + c) / len(sequence) * 100
    return gc_content

def gc_content_of_fasta(fasta_file):
    """
    Calculate the GC content for each sequence in a FASTA file.
    
    Parameters:
    fasta_file (str): Path to the FASTA file.
    
    Returns:
    list: A list of GC content values for each sequence.
    """
    gc_contents = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        gc_contents.append(calculate_gc_content(str(record.seq)))
    return gc_contents

def process_fasta_files_in_directory(directory):
    """
    Process all FASTA files in a directory to calculate average GC content.
    
    Parameters:
    directory (str): Path to the directory containing FASTA files.
    
    Returns:
    list: A list of tuples containing filename and average GC content.
    """
    results = []
    for filename in os.listdir(directory):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            filepath = os.path.join(directory, filename)
            gc_contents = gc_content_of_fasta(filepath)
            if gc_contents:
                average_gc_content = sum(gc_contents) / len(gc_contents)
                results.append((filename, average_gc_content))
    return results

def save_results_to_csv(results, output_file):
    """
    Save the GC content results to a CSV file.
    
    Parameters:
    results (list): A list of tuples containing filename and average GC content.
    output_file (str): Path to the output CSV file.
    """
    with open(output_file, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Filename", "Average GC Content (%)"])
        for row in results:
            writer.writerow(row)

# Directory containing FASTA files
directory = '.'

# Output CSV file
output_file = 'average_gc_content_results.csv'

# Process FASTA files and save results to CSV
results = process_fasta_files_in_directory(directory)
save_results_to_csv(results, output_file)

print(f"Average GC content results have been saved to {output_file}")
