import sys
from Bio.Seq import Seq

def load_sequence(file_path):
    """
    Reads the DNA sequence from a text file.
    Removes whitespace and newlines.
    """
    try:
        with open(file_path, "r") as file:
            sequence = file.read().replace("\n", "").strip().upper()
        return sequence
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        sys.exit(1)

def dna_to_protein(dna_sequence):
    """
    Translates a DNA sequence into a protein sequence.
    """
    seq = Seq(dna_sequence)
    # Translate DNA to Protein (to_stop=False keeps stop codons as *)
    protein_seq = seq.translate()
    return str(protein_seq)

def generate_variants(dna_sequence):
    """
    Generates all possible single nucleotide variants.
    Filters out synonymous mutations and stop codons.
    """
    original_protein = dna_to_protein(dna_sequence)
    variants = []
    seen_proteins = set()
    
    bases = ['A', 'T', 'C', 'G']
    stop_codon_count = 0
    synonymous_count = 0
    
    print("Starting simulation...")
    
    for i in range(len(dna_sequence)):
        original_base = dna_sequence[i]
        
        for base in bases:
            if base == original_base:
                continue
            
            # Create mutation
            mutated_dna = dna_sequence[:i] + base + dna_sequence[i+1:]
            mutated_protein = dna_to_protein(mutated_dna)
            
            # Filter 1: Check for Stop Codons (*)
            if "*" in mutated_protein[:-1]: # Ignore stop at the very end
                stop_codon_count += 1
                continue
            
            # Filter 2: Check for Synonymous Mutations (No amino acid change)
            if mutated_protein == original_protein:
                synonymous_count += 1
                continue
            
            # Filter 3: Check for Duplicates
            if mutated_protein not in seen_proteins:
                variants.append({
                    'position': i + 1,
                    'original_base': original_base,
                    'new_base': base,
                    'protein_seq': mutated_protein
                })
                seen_proteins.add(mutated_protein)

    return original_protein, variants, stop_codon_count, synonymous_count

def save_results(original_protein, variants, output_file="mutation_results.txt"):
    """
    Writes the comparison of variants to a text file.
    Uses ^ to mark changes.
    """
    with open(output_file, "w") as outfile:
        outfile.write(f"Original Sequence:\n{original_protein}\n\n")
        outfile.write("-" * 50 + "\n")
        outfile.write(f"Total Unique Variants Found: {len(variants)}\n")
        outfile.write("-" * 50 + "\n\n")
        
        for idx, variant in enumerate(variants):
            outfile.write(f"Variant {idx+1} (Nt Position: {variant['position']} | {variant['original_base']}->{variant['new_base']}):\n")
            
            # Visual comparison logic
            mutated_seq = variant['protein_seq']
            display_seq = ""
            for j in range(len(original_protein)):
                if original_protein[j] == mutated_seq[j]:
                    display_seq += mutated_seq[j]
                else:
                    # Mark the change with special characters
                    display_seq += f"[{mutated_seq[j]}]" 
            
            outfile.write(f"{display_seq}\n\n")
            
    print(f"Results successfully saved to {output_file}")

def main():
    # File name containing the wild type sequence
    input_file = "input_sequence.txt"
    
    # 1. Load Data
    dna_seq = load_sequence(input_file)
    print(f"Loaded DNA Sequence Length: {len(dna_seq)} bp")
    
    # 2. Process Variants
    wt_protein, all_variants, stops, silent = generate_variants(dna_seq)
    
    # 3. Print Statistics (Reference for your thesis)
    print("\n--- Summary Statistics ---")
    print(f"Wild Type Protein Length: {len(wt_protein)} aa")
    print(f"Total Theoretical Variants: {len(all_variants)}")
    print(f"Stop Codons Eliminated: {stops}")
    print(f"Synonymous Mutations Eliminated: {silent}")
    
    # 4. Save to File
    save_results(wt_protein, all_variants)

if __name__ == "__main__":
    main()
