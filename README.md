# DNA Palindrome Explorer

## Overview
DNA Palindrome Explorer is an interactive web application built using Streamlit. It analyzes a given DNA sequence to identify palindromic regions by distinguishing between two types:

- **Split Palindromes:** Occur as two separate segments in the DNA that are reverse complements of each other.
- **Direct Palindromes:** Continuous sequences that are identical to their reverse complement.

DNA palindromes are significant in understanding gene regulation, DNA structure, and evolutionary processes.

## Features
- **Multiple Input Options:**  
  - **Manual Input:** Directly type your DNA sequence.
  - **FASTA File Upload:** Upload a FASTA file; if multiple sequences are detected, choose the one you want to analyze.
- **Customizable Analysis Parameters:**  
  - Set the minimum palindrome length to tailor the analysis.
- **Backend Logic:**
  - **Input Handling:** Processes both manually entered sequences and FASTA file inputs.
  - **Sequence Validation:** Ensures the input contains only valid nucleotide bases (A, T, C, G).
  - **Reverse Complement Calculation:** Converts the sequence into its reverse complement based on standard base pairing rules.
  - **Palindrome Detection:**
    - Scans the sequence for all possible subsequences (of at least a given length) that are palindromic.
    - Filters and categorizes the results into **Split Palindromes** and **Direct Palindromes**.
  - **Output Formatting:**  
    - Displays the input sequence in a user-friendly format with 60 bp lines (segmented in 10 bp groups).
    - Summarizes the detected palindromes with their lengths, positions, and counts.
