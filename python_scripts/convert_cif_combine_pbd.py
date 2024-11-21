import os
from Bio.PDB import MMCIFParser, PDBIO, Chain

# Set the directory where your .cif files are located
cif_directory = './'  # Update this path to where your files are
output_pdb = 'combined_structure.pdb'  # Name of the output PDB file

# Initialize the parser and the PDBIO writer
parser = MMCIFParser(QUIET=True)
io = PDBIO()

# Initialize an empty structure
combined_structure = None

# Loop through all files in the directory
for file_name in sorted(os.listdir(cif_directory)):
    if file_name.endswith(".cif"):
        cif_file_path = os.path.join(cif_directory, file_name)
        
        # Parse the CIF file
        structure = parser.get_structure(file_name, cif_file_path)
        
        # Combine the structure into the main structure
        if combined_structure is None:
            combined_structure = structure
        else:
            # Add chains from the new structure to the combined structure
            for model in structure:
                for chain in model:
                    # Append a unique chain ID if needed
                    chain.id = chr(ord(chain.id) + len(combined_structure[0])) # to avoid chain ID conflicts
                    combined_structure[0].add(chain)

# Save the combined structure as a PDB file
io.set_structure(combined_structure)
io.save(os.path.join(cif_directory, output_pdb))

print(f"Combined PDB file saved as {output_pdb}")
