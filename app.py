import streamlit as st
import io

# Set page configuration to wide layout and define a simple page title.
st.set_page_config(layout="wide", page_title="DNA Palindrome Explorer")


# ----------------------
# Utility Functions
# ----------------------

def manual_input_sequence(input_text: str) -> str:
    """Process manual text input."""
    return input_text.strip().lower()


def fasta_input(file_content: bytes) -> dict:
    """Parse FASTA file content and return a dictionary mapping headers to sequences."""
    content = file_content.decode("utf-8")
    lines = content.splitlines()
    sequence_header = []
    sequence_lines = []
    sub_lines = []
    for idx, line in enumerate(lines):
        line = line.strip()
        if line.startswith(">"):
            if sub_lines:
                sequence_lines.append(''.join(sub_lines).lower())
                sub_lines = []
            sequence_header.append(line)
        else:
            sub_lines.append(line)
            if idx == len(lines) - 1 and sub_lines:
                sequence_lines.append(''.join(sub_lines).lower())
    sequences = dict(zip(sequence_header, sequence_lines))
    return sequences


def valid_sequence(sequence: str) -> bool:
    valid_base = {'a', 't', 'c', 'g'}
    return all(base in valid_base for base in sequence)


def reverseComp(sequence: str) -> str:
    rules_for_complement = {'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    complement_seq = "".join(rules_for_complement.get(base, base) for base in sequence)
    return complement_seq[::-1]


def all_possible_palindromes(sequence: str, reverseComp_seq: str, minLength: int):
    seqLength = len(sequence)
    palindromeList = []
    for length in range(minLength, seqLength):
        for start_idx in range(seqLength - length + 1):
            sub_sequence = sequence[start_idx:start_idx + length]
            if sub_sequence in reverseComp_seq and sub_sequence not in [p[1] for p in palindromeList]:
                palindromeList.append((start_idx + 1, sub_sequence, start_idx + length))
    return palindromeList if palindromeList else None


def find_longest(all_pals):
    longest_palindrome = []
    for i in range(len(all_pals)):
        current_sub_sequence = all_pals[i]
        current_start, current_end = current_sub_sequence[0], current_sub_sequence[2]
        flag = True
        for j in range(len(all_pals)):
            if i != j:
                compare_sub_sequence = all_pals[j]
                compare_start, compare_end = compare_sub_sequence[0], compare_sub_sequence[2]
                if (current_start >= compare_start) and (current_end <= compare_end):
                    flag = False
                    break
        if flag:
            longest_palindrome.append(current_sub_sequence)
    return longest_palindrome


# Renamed "spacer" to "split_palindromes" for clarity.
def spacer(longest_pals):
    split_dict = {}
    for item in longest_pals:
        current_start, current_segment, current_end = item
        rev_current = reverseComp(current_segment)
        for compare_item in longest_pals:
            compare_start, compare_segment, compare_end = compare_item
            if current_start != compare_start and rev_current == compare_segment and item not in split_dict.values():
                split_dict[item] = compare_item
    return split_dict


# Revised nonspacer function: Direct palindromes are those equal to their full reverse complement.
def nonspacer(all_pals):
    direct_list = []
    for item in all_pals:
        _, sub_sequence, _ = item
        if sub_sequence == reverseComp(sub_sequence):
            direct_list.append(item)
    return direct_list


def format_output(sequence: str, min_length: int, split_dict, direct_list) -> str:
    seqLength = len(sequence)
    # Format the query sequence in 60 bp lines with 10 bp sub-chunks.
    query_list = []
    query_bp_count = 1
    for i in range(0, seqLength, 60):
        query_line = sequence[i:i + 60]
        segmented_line = " ".join(query_line[j:j + 10] for j in range(0, len(query_line), 10))
        query_list.append(f'{query_bp_count:>4} {segmented_line}')
        query_bp_count += 60
    query_string = "\n".join(query_list)

    # Format split (formerly spacer) palindrome display.
    split_output = ["-"] * seqLength
    for (start1, seg1, end1), (start2, seg2, end2) in split_dict.items():
        split_output[start1 - 1:end1] = list(seg1)
        split_output[start2 - 1:end2] = list(seg2)
    split_list = []
    split_bp_count = 1
    for i in range(0, seqLength, 60):
        line = split_output[i:i + 60]
        segmented_line = " ".join("".join(line[j:j + 10]) for j in range(0, len(line), 10))
        split_list.append(f'{split_bp_count:>4} {segmented_line}')
        split_bp_count += 60
    split_string = "\n".join(split_list)

    # Format direct (formerly non-spacer) palindrome display.
    direct_output = ["-"] * seqLength
    for (start, seg, end) in direct_list:
        direct_output[start - 1:end] = list(seg)
    direct_list_lines = []
    direct_bp_count = 1
    for i in range(0, seqLength, 60):
        line = direct_output[i:i + 60]
        segmented_line = " ".join("".join(line[j:j + 10]) for j in range(0, len(line), 10))
        direct_list_lines.append(f'{direct_bp_count:>4} {segmented_line}')
        direct_bp_count += 60
    direct_string = "\n".join(direct_list_lines)

    output = []
    output.append("DNA Palindrome Explorer Results\n\n")
    output.append(f"Query sequence ({seqLength} bp):\n{query_string}\n\n")
    output.append(f"Analysis for Palindromes longer or equal to {min_length} bp\n\n")

    if split_dict:
        output.append(f"Split Palindrome(s), count({len(split_dict)}):\n")
        output.append(split_string + "\n\n")
        output.append("Split palindrome lengths & positions:\n")
        for (start1, seg1, end1), (start2, seg2, end2) in split_dict.items():
            output.append(f'{len(seg1):>2} bp -> [{start1} {seg1} {end1}] : [{start2} {seg2} {end2}]\n')
        output.append("\n")
    else:
        output.append("No split palindromes detected.\n\n")

    if direct_list:
        output.append(f"Direct Palindrome(s), count({len(direct_list)}):\n")
        output.append(direct_string + "\n\n")
        output.append("Direct palindrome lengths & positions:\n")
        for start, seg, end in direct_list:
            output.append(f'{len(seg):>2} bp -> [{start} {seg} {end}]\n')
    else:
        output.append("No direct palindromes detected.\n")

    return "".join(output)


# ----------------------
# Streamlit App Layout
# ----------------------

# Sidebar for user inputs to keep the main area clean.
with st.sidebar:
    st.header("Input Settings")
    input_method = st.radio("Select input method", ("Manual input", "Upload FASTA file"))
    min_length = st.number_input("Minimum palindrome length", min_value=2, value=4, step=1)

    # Manual input text area in sidebar if chosen.
    if input_method == "Manual input":
        manual_text = st.text_area("Enter the gene sequence:", height=150)
    else:
        uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "txt"])

# Main section with explanation and results.
st.markdown("""
# DNA Palindrome Explorer

**What are DNA Palindromes?**  
In genetics, a DNA palindrome is a segment of DNA in which the sequence of nucleotides is identical when read from 5' to 3' on one strand and from 5' to 3' on the complementary strand. This means that the sequence on one strand is the reverse complement of the sequence on the other strand. DNA palindromes are not true palindromes in the linguistic sense because the sequence of bases on one strand is not exactly the same as the sequence on the other strand. Instead, they are reverse complements.

For example, consider the DNA sequence 5'-CATATC-3'. Its complementary strand would be 3'-GTATAC-5'. If you read the complementary strand from 5' to 3', it becomes 5'-CATATC-3', which is the same as the original strand. Therefore, CATATC is a DNA palindrome. These palindromic sequences are often the recognition sites of restriction enzymes, proteins that can cut DNA at specific sequences.

**Types of Palindromes:**
- **Split Palindromes:** Occur as two separate segments that are reverse complements of each other.
- **Direct Palindromes:** Continuous sequences that are identical to their reverse complement.
""")

sequence = None
if input_method == "Manual input":
    if 'manual_text' in locals() and manual_text:
        sequence = manual_input_sequence(manual_text)
else:
    if 'uploaded_file' in locals() and uploaded_file is not None:
        try:
            content = uploaded_file.read()
            sequences = fasta_input(content)
            if len(sequences) == 1:
                sequence = list(sequences.values())[0]
                st.info("Single sequence detected from file.")
            elif len(sequences) > 1:
                seq_options = {header: seq for header, seq in sequences.items()}
                selected_header = st.selectbox("Multiple sequences detected. Select one:", list(seq_options.keys()))
                sequence = seq_options[selected_header]
        except Exception as e:
            st.error("Error reading FASTA file: " + str(e))

if st.button("Run Analysis"):
    if not sequence:
        st.error("No sequence provided.")
    elif not valid_sequence(sequence):
        st.error("Invalid sequence. Ensure it contains only A, T, C, G characters.")
    else:
        reverse_sequence = reverseComp(sequence)
        all_pals = all_possible_palindromes(sequence, reverse_sequence, min_length)
        if all_pals:
            longest_pals = find_longest(all_pals)
            # Using updated naming: "split" for split palindromes, and "direct" for direct palindromes.
            split_palindromes = spacer(longest_pals)
            direct_palindromes = find_longest(nonspacer(all_pals))
            output_str = format_output(sequence, min_length, split_palindromes, direct_palindromes)
            st.code(output_str, language="plaintext")
        else:
            st.info("No palindromes detected in the provided sequence.")