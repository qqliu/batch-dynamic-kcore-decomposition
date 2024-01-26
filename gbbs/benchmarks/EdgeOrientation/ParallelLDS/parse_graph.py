# Import sys module to access command line arguments
import sys

# Check if the user provided two file names as arguments
if len(sys.argv) != 4:
    print("Usage: python program.py input_file output_file")
    exit()

# Assign the input and output file names to variables
input_file = sys.argv[1]
output_file = sys.argv[2]
output_file2 = sys.argv[3]

# Create a set to store the seen words
seen_words = set()
counter = 0
vertices = {}

# Open the input file in read mode and the output file in write mode
with open(input_file, "r") as fin, open(output_file, "w") as fout, open(output_file2, "w") as fout2:
    # Loop through each line in the input file
    for line in fin:
        # Split the line into a list of words
        words = line.split()
        # Remove the first and last words from the list
        modified_words = words[1:-1]
        # Join the modified words with spaces
        modified_line = " ".join(modified_words)
        # Check if the modified line has already been seen
        new_edge = (int(modified_words[0]), int(modified_words[1]))
        if new_edge not in seen_words:
            u = int(modified_words[0])
            if u not in vertices:
                vertices[u] = counter
                counter += 1
            v = int(modified_words[1])
            if v not in vertices:
                vertices[v] = counter
                counter += 1

            edge = (u, v)
            edge_twin = (v, u)

            modified_line = str(vertices[u]) + " " + str(vertices[v])
            # Write the modified line to the output file
            fout.write(modified_line + "\n")
            fout2.write("+ " + modified_line + " 1\n")
            # Add the modified line to the seen words set
            seen_words.add(edge)
            seen_words.add(edge_twin)
        else:
            print("DUPLICATE FOUND")

