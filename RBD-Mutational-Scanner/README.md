# SARS-CoV-2 RBD Deep Mutational Scanning Analysis

This repository contains the bioinformatic scripts used to determine the theoretical mutational landscape of the SARS-CoV-2 Spike protein Receptor Binding Domain (RBD).

This analysis is part of the PhD research conducted at the University of Leicester.

## Overview

The script `rbd_variant_scanner.py` performs a comprehensive computational scan of a given DNA sequence. It simulates every possible single nucleotide substitution and translates the resulting sequences to proteins.

### Key Features:
* **Translation:** Converts DNA to Protein using Biopython.
* **Mutation Simulation:** Iterates through every nucleotide position (A, T, G, C).
* **Filtering:** Automatically excludes:
    * Synonymous (silent) mutations.
    * Nonsense mutations (Stop codons).
    * Duplicate protein sequences.
* **Output:** Generates a detailed report highlighting amino acid changes.

## Requirements

* Python 3.x
* Biopython

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
   ```

2. Install dependencies:
   ```bash
   pip install biopython
   ```

## Usage

1. Place your target DNA sequence in `input_sequence.txt`.
2. Run the script:
   ```bash
   python rbd_variant_scanner.py
   ```
3. Check the output in `mutation_results.txt` and the summary statistics in the console.

## Results Reference

As cited in the thesis, this tool identified **1,323 distinct alternative protein sequences** for the RBD gene, while filtering out **75 stop-codon mutations**.

## License

[MIT License](https://opensource.org/licenses/MIT)
