'''
The Python script below analyses gene palindrome(s) exist in the base pairs of a given query sequence,
whether in the form of spatially separated or non-spatially separated palindrome.

Pseudocodes:
1. Request sequence input (manual or FASTA format) from user.
    1.1 If FASTA input has multiple sequence, request user which sequence to perform analysis on.
2. Check if query sequence is a valid DNA sequence.
3. Request minimum palindrome length of interest from user.
4. Convert query sequence to its reverse-complement counterpart.
5. Compare the query sequence with reverse-complemented sequence starting with a defined minimum length,
   then increment the length by 1 bp for each iteration until maximum length of sequence.
   This will return all possible palindromes exist in the DNA sequence, ranging from the defined minimum length
   to the maximum length of palindrome.
6. From the list of all palindromes identified in the DNA sequence,
    identify spacer palindrome and non-spacer palindromes:-
    6.1. Spacer palindrome:
        From the list of all palindromes found in sequence, isolate the longest palindromic sequence from their
        short, repeated palindromic sequence. Then, identify spacer palindrome by means of a) the difference of start
        (or end) index number; and b) the reverse-complement counterparts separated spatially in the sequence.
        The matched palindrome output will return as a key-value pair in the form of dictionary.
    6.2. Non-spacer palindrome:
        From the list of all palindromes found in sequence, identify all non-spacer palindrome by means of a) their
        start (or end) index numbers must be equal; b) the palindromic sequence is the same when subject to reverse-
        complement conversion of bases. Then, filter the short, repeated palindromes to avoid redundancy.
        The palindrome output will return as item in list.
7. Generate an output file to present analysis result.
'''

# Function to extract sequence from manual input
def manual_input():
    print("Enter the gene sequence: ")
    sequence = input("> ")
    sequence = sequence.lower()

    return sequence

# Function to extract sequence from FASTA file format
# This function assumes FASTA file may contain more than one DNA sequences,
# it prompts user to choose only one sequence to proceed with analysis
def fasta_input():
    print("Enter FASTA file name:")
    fasta_file = input("> ")

    # Open fasta file based on user input
    with open(fasta_file, 'r') as file:
        lines = file.readlines() # Separate each line of fasta file into list
        sequence_header = [] # Initialize list to store header of sequences
        sequence_lines = [] # Initialize list to store different sequences
        sub_lines = [] # Initialize list to temporarily store one sequence

        # Individually check each line of lines
        for line in lines:
            # Conditional statement for last line of fasta file; store last seq
            if not line.startswith(">") and line == lines[len(lines)-1]:
                sub_lines.append(line.strip())
                sub_seq = ''.join(sub_lines).lower()
                sequence_lines.append(sub_seq)

            # Conditional statement for line NOT header, append line into sub_lines
            elif not line.startswith(">"):
                sub_lines.append(line.strip())

            # Conditional statement for new sequence header; store previous seq into sequence_lines
            elif line.startswith(">") and sub_lines != []:
                sub_seq = ''.join(sub_lines).lower()
                sequence_lines.append(sub_seq)
                sub_lines = []

            # Conditional statement for header, store header into sequence_header
            if line.startswith(">"):
                sequence_header.append(line.strip())

    # Create dictionary using sequence_header as KEYS and sequence_lines as VALUES
    sequences = dict(zip(sequence_header,sequence_lines))

    if sequences:
        return sequences

    else:
        print("Error reading FASTA file: Please try again.")
        return None
        # Returns: The gene sequence extracted from the FASTA file, or an error message if the file is invalid.

# Checks if the sequence contains only valid nucleotide bases.
# This validation is crucial to avoid errors during the palindrome analysis.
def valid_sequence(sequence):
    valid_base = ['a', 't', 'c', 'g']

    for base in sequence:
      if base not in valid_base:
        return None

    return sequence

def user_input():
    # Function in user_input() to RETURN valid sequence
    def return_valid(sequence):
        if valid_sequence(sequence) is None:
            print("Invalid input. The sequence should only contain 'A', 'T', 'G', and 'C'.")
            return user_input()
        else:
            return sequence

    prompt = input("Please enter 'm' to input sequence manually, 'f' to use FASTA format, or 'q' to quit (m/f/q):\n> ")
    # This guides the user to input the sequence either manually or via a FASTA file.

    if prompt.lower() == 'm':
        sequence = manual_input()
        return return_valid(sequence) # Return if sequence is valid

    elif prompt.lower() == 'f':
        sequences = fasta_input()
        keys = list(sequences.keys())

        if len(sequences.keys()) == 1: # Check dictionary if sequences has only one sequence
            sequence = sequences[keys[0]] # Assign the single line into sequence
            return return_valid(sequence) # Return valid sequence
        else:
            print ('Multiple sequences have been detected in file.')
            i = 0 # Counter for index
            for key in sequences:
              print (f'({i+1}): {key}, {sequences[key][:10]}...')
              i += 1

            # While loop to keep looping if user does not input index properly
            while (i != 0):
                # Store user input index for sequence
                seqIndex = input("Type the number of which sequence to perform Palindrome analysis:\n> ")
                if seqIndex.isdigit(): # Check if user input is a number
                    i = 0 # Change flag value to exit loop
                    sequence = sequences[keys[int(seqIndex)-1]] # assign sequence of user's index
                else:
                  print("ERROR: Please input a number of the sequence.\n") # Error statement

            return return_valid(sequence) #return value

    elif prompt.lower() == 'q':
        return None

    else:
        print("Invalid input. Please enter 'm', 'f', or 'q'.")
        return user_input()
        # Incorporates a recursive call to handle invalid inputs and ensure correct data entry.

# Reverse-complement conversion of bases in DNA sequence
def reverseComp(sequence):
    rules_for_complement = {'a':'t', 't':'a', 'g':'c', 'c':'g'} # It adheres to the DNA base pairing rules: A-T and G-C.
    complement_seq = "" # Initiate an empty string

    for i in range(0, len(sequence)): # Iterate in the length of sequence
        base = sequence[i] # Variable to store iterative base
        complement_seq += rules_for_complement[base] # Obtain rules items then add to complement_seq string
    reverseComp_seq = complement_seq[::-1] # Reverse the complement_seq string

    return reverseComp_seq # Return the reverse complement of the input sequence

# Find all possible palindrome in a query sequence
# This function is crucial for identifying regions in DNA that are symmetric and may have biological significance.
def all_possible_palindromes(sequence, reverseComp_seq, minLength):
    # Input requirements:
    # - sequence: The original DNA sequence.
    # - reverseComp_seq: The reverse complement of the DNA sequence, which is used to identify palindromes.
    # - minLength: The minimum length of palindrome sequences to be considered.

    seqLength = len(sequence) # Compute the sequence length
    palindromeList = [] # Initialize an empty list to store the details of identified palindrome sequences.

    # Define how many base pair check per iteration through specified minLength,
    # then loop through the entire length of sequence
    for length in range(minLength, seqLength - 1):

        # Iterates over all possible starting indices for subsequences of a specified length within a range of given sequence
        # Deduction of length is needed to prevent loop goes beyond boundaries
        for start_idx in range(seqLength - length + 1):
            # For each loop, extracts a subsequence from the original sequence
            sub_sequence = sequence[start_idx:start_idx+length]

            # If a palindrome is found and is not already in the list, add its details (start index, sequence, end index).
            if sub_sequence in reverseComp_seq and sub_sequence not in palindromeList:
                palindromeList.append((start_idx+1, sub_sequence, start_idx+len(sub_sequence)))

    minLength += 1 # Increase the minimum length by 1 bp for the next iteration.

    if palindromeList:
        return palindromeList
    # A list of tuples, each containing the start index, palindrome sequence, and end index of each identified palindrome.

    else:
        return None # If no palindrome is found, return None.


# This function is important for distinguishing the most significant palindrome sequences by utilizing index numbers,
# as shorter ones might be part of longer palindromes.
def find_longest(all_pals):
# all_pals: A list of tuples containing details of all identified palindrome sequences (start index, sequence, end index).
    # First initialize an empty list to store the longest unique palindrome sequences.
    longest_palindrome = []

    for i in range(len(all_pals)):
        current_sub_sequence = all_pals[i]
        # For each sequence in current_sub_sequence, extract its start and end indexes.
        current_start, current_end = current_sub_sequence[0], current_sub_sequence[2]

        # Initiate a flag to control the conditional statement, also assume every sequence is unique in the first place.
        flag = True
        for j in range(len(all_pals)):
            if i != j: # If current_sub_sequence index is not compare_sub_sequence index
                # Start comparing every other sequence in the list except for current_sub_sequence.
                compare_sub_sequence = all_pals[j]
                # For each sequence in compare_sub_sequence, extract its start and end indexes
                compare_start, compare_end = compare_sub_sequence[0], compare_sub_sequence[2]

                if (current_start <= compare_start) and (current_end >= compare_end):
                    continue # If the current sequence entirely encompasses another, flag remains true.

                elif (current_start >= compare_start) and (current_end <= compare_end):
                    flag = False # If the current sequence is smaller than the other, flag the sequence as false.
                    break

        if flag:
            # Append all the true flagged sequence as the longest and/or unique palindrome(s)
            longest_palindrome.append(current_sub_sequence)

    return longest_palindrome

# Identify spacer palindromes and map the palindrome sequence with its counterpart
def spacer(longest_pals):
# longest_pals: This function requires the longest/unique palindromes identified from list of all palindromes.
    spacer = {} # Initialize an empty dictionary to store spacer palindromes.

    for item in longest_pals:
        # For each palindrome (item), extract its start index, sequence, and end index.
        current_start, current_segment, current_end = item
        reverseComp_segment = reverseComp(current_segment) # Generate the reverse complement of the palindrome sequence.

        # For compare item (the rest of palindromes in the list), extract its start index, sequence, and end index.
        for compare_item in longest_pals:
            compare_start, compare_segment, compare_end = compare_item

            # If the current palindrome and another palindrome are located at different positions
            # and are reverse complements of each other, and the current palindrome is not already included
            # in the spacer dictionary as value (to avoid double counting), consider it as a spacer palindrome.
            if current_start != compare_start and reverseComp_segment == compare_segment and \
            item not in spacer.values():
                # Add the current palindrome and its matching palindrome as a key-value pair in the spacer dictionary.
                spacer[item] = compare_item

    return spacer
    # A dictionary of spacer palindromes. Each key-value pair represents a spacer palindrome and its counterpart in the sequence.

# Find nonspacer from list of all possible palindromes
def nonspacer(all_pals):
# all_pals: A list of tuples containing details of all identified palindrome sequences (start index, sequence, end index).
    nonspacer = []  # Initialize an empty list to store non-spacer palindromes.

    for i in range(len(all_pals)):
        current_sub_sequence = all_pals[i] # Iterate through each palindrome item in the provided list.
        sub_sequence = current_sub_sequence[1] # For each item, extract its palindrome sequence (located at index 1).

        complement_sub_sequence = "" # Initiate an empty string to store converted sequence
        for base in reversed(sub_sequence):   # Iterate over the bases of the current segment in reverse order.
            # For each base, find its complement and add it to a string.
            complement_sub_sequence += reverseComp(base)

        complement_sub_sequence = "".join(complement_sub_sequence)

        # If they are equal, it indicates the sequence is a non-spacer palindrome (self-complementary within the same region).
        if sub_sequence == complement_sub_sequence:
            nonspacer.append(current_sub_sequence)  # Add the non-spacer palindrome to the list.

    return nonspacer  # Return the list containing non-spacer palindromes.

# Generate an external file to present analysis result
# This function is for documenting and visualizing the findings in a clear, readable format.
def output_file(filename, sequence, min_length, spacer, nonspacer):
    # 1. Initialize graphical representations for the sequence, spacer, and non-spacer palindromes.
    seqLength = len(sequence)
    spacer_output = ["-"] * seqLength # Initializes a list with multiple(seqLength) '-' characters.
    nonspacer_output = ["-"] * seqLength # Similar to spacer_output, but for non-spacer palindromes.

    # 2. Construct a formatted string representation of the query sequence, segmenting it for readability.
    query_bp_count = 1 # A counter for base pair positions in the sequence, starting from 1.
    query_list = [] # An empty list to hold formatted lines of the sequence.

    # Extract a 60 base pair long segment from the entire length of sequence
    for i in range(0, seqLength, 60):
        query_line = sequence[i:i + 60]

        # Further divide this segment into 10 base pair chunks for readability, separating them with spaces.
        segmented_line = ""
        for j in range(0, len(query_line), 10):
            segmented_line += query_line[j:j + 10] + " "

        # Append the formatted line, along with the starting base pair position, to the query_list.
        query_list.append(f'{query_bp_count:>4} {segmented_line}')
        query_bp_count += 60

    # Joins all formatted lines from query_list into a single string, separated by new lines.
    query_string = '\n'.join(query_list)

    # 3. Create a graphical representation of spacer palindromes
    #   a. Mark spacer palindrome positions within the sequence.
    #   b. Format and align these representations with the query sequence.
    spacer_bp_count = 1
    spacer_list = []

    # For each spacer palindrome, mark its sequence on the 'spacer_output' list at the corresponding positions
    for (start1, spacer_segment1, end1), (start2, spacer_segment2, end2) in spacer.items():

        # This involves replacing '-' characters with the actual bases of the spacer palindrome sequence.
        spacer_output[start1 - 1:end1] = list(spacer_segment1)
        spacer_output[start2 - 1:end2] = list(spacer_segment2)
        # In Py index starts for 0. But in seq index starts from 1. So it's needed to "-1"

    # Segment the 'spacer_output' list into 60 base pair chunks, further divided into 10 base pair sub-chunks.
    for i in range(0, seqLength, 60):
        spacer_line = spacer_output[i:i + 60]
        segmented_line = []

        for j in range(0, len(spacer_line), 10):
            segmented_line.append(''.join(spacer_line[j:j + 10]) + ' ')

        # Append each formatted line to 'spacer_list' with the corresponding base pair position.
        segmented_line = ''.join(segmented_line)
        spacer_list.append(f'{spacer_bp_count:>4} {segmented_line}')
        spacer_bp_count += 60

    spacer_string = '\n'.join(spacer_list) # Combine two different format into a single string for output

    # 4. Create a graphical representation of non-spacer palindromes:
    #   a. Mark non-spacer palindrome positions within the sequence.
    #   b. Format and align these representations with the query sequence.
    #   Similar method as shown above
    nonspacer_bp_count = 1
    nonspacer_list = []

    for (start, nonspacer_segment, end) in nonspacer:
        nonspacer_output[start - 1:end] = list(nonspacer_segment)

    for i in range(0, seqLength, 60):
        nonspacer_line = nonspacer_output[i:i + 60]
        segmented_line = []

        for j in range(0, len(nonspacer_line), 10):
            segmented_line.append(''.join(nonspacer_line[j:j + 10]) + ' ')

        segmented_line = ''.join(segmented_line)
        nonspacer_list.append(f'{nonspacer_bp_count:>4} {segmented_line}')
        nonspacer_bp_count += 60

    nonspacer_string = '\n'.join(nonspacer_list)

    # 5. Write the formatted analysis results to the specified file:
    #   a. Include the query sequence and its length.
    #   b. Detail the analysis of palindromes, including counts, palindrome lengths and positions.
    #   c. Present spacer and non-spacer palindrome information separately for clarity.
    with open(filename, 'w') as f:
        f.write("Palindromes Detection in DNA Sequence" + "\n" * 2) # File's header

        # Include the formatted query sequence with its length.
        f.write(f"Query sequence ({seqLength} bp):\n")
        f.write(query_string)
        f.write("\n" * 2)

        # Informing user the specified minimum length for palindromes
        f.write(f"Analysis for Palindromes longer or equal to {min_length} bp" + "\n" * 2)

        if spacer: # If spacer palindrome(s) detected
            f.write(f"Spacer Palindrome(s), count({len(spacer)}):\n") # Compute the total count of spacer palindromes
            # Include the graphical representation of spacer palindromes aligned with the sequence.
            f.write(spacer_string + '\n' * 2)
            # List the length and positions of each spacer palindrome for detailed analysis.
            f.write("Spacer palindrome length & positions:\n")
            for (start1, spacer_segment1, end1), (start2, spacer_segment2, end2) in spacer.items():
                f.write(f'{len(spacer_segment1):>2} bp -> [{start1} {spacer_segment1} {end1}] : [{start2} {spacer_segment2} {end2}]\n')
            f.write("\n")
        else:
            f.write("Spacer palindrome not detected." + "\n"*2) # A statement indicates that no spacer palindromes were detected.

        if nonspacer: # Similar method as shown above
            f.write(f"Non-spacer Palindrome(s), count({len(nonspacer)}):\n")
            f.write(nonspacer_string + '\n' * 2)
            f.write("Non-spacer palindrome length & positions:\n")
            for start, nonspacer_segment, end in nonspacer:
                f.write(f'{len(nonspacer_segment):>2} bp -> [{start} {nonspacer_segment} {end}]\n')
        else:
            f.write("Non-spacer palindrome not detected.")

# Main function to assemble all the functions and execute them in different steps.
def main():
    print("Welcome to the DNA Palindrome Analyzer Tool!\n")
    # Requesting sequence input from user
    sequence = user_input()

    if sequence is not None: # Check if a valid sequence is provided.
        # Requesting minimum/initial length of palindrome needed for analysis from user
        min_length = int(input("Please enter the minimum length of the palindrome you are interested in: \n> "))

        reverseComp_sequence = reverseComp(sequence) # Generate the reverse-complement sequence
        # Identify all possible palindromes in the sequence
        all_pals = all_possible_palindromes(sequence, reverseComp_sequence, min_length)

        # If all_possible_palindrome() does not return None, execute the following functions
        if all_pals:
            # Spacer: First find the longest, unique palindromes in all_pals list, then identify spacer palindromes
            # 1. Find the longest unique palindromes
            longest_pals = find_longest(all_pals)
            # 2. Identified spacer palindromes and map to their reverse-complement counterparts
            spacer_pals = spacer(longest_pals)

            # Non-spacer: First find all the non-spacer palindrome in all_pals list,
            # then isolate the longest palindrome from short redundancy palindromes
            # 1. Find all the non-spacer palindromes
            nonspacer_pals = nonspacer(all_pals)
            # 2. Filter out the longest non-spacer palindromes from their shorter versions
            nonspacer_pals = find_longest(nonspacer_pals)

            # Generate an output file
            filename = 'PalindromeAnalysis.txt'
            # Write all details into external file
            output_file(filename, sequence, min_length, spacer_pals, nonspacer_pals)

            # Following lines provide a summary of the findings to the user and print it to terminal
            print("")
            print("Palindrome(s) found in the query sequence.")
            print(f"Spacer palindrome count: ({len(spacer_pals)})")
            print(f"Non-spacer palindrome count: ({len(nonspacer_pals)})")
            print(f"Analysis file '{filename}' has been generated.")

        else:
            # Handle the case where no palindromes are detected.
            print("")
            print("No palindrome is detected in the query sequence.")

if __name__ == "__main__":
    main()