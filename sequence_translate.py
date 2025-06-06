######################## BCB 546X: Python Assignment Details ########################

# ** Your Mission: Complete Python code in a Jupyter Notebook ** #

# In these commented lines, you will find a list of tasks. Below, within the code, you will find
# suggestions on how to complete this. Conveniently, you can load this script into a Jupyter cell
# using the `%load` magic command. You might find that useful. 

#-- Functions --#
## 1. Document Dr. X's function with comments and with markdown text in your Jupyter notebook.
## 2. Write a function that translates a string of nucleotides to amino acids based on Dr. X's pseudo-code suggestion.
## 3. Write an alternative translation function.
## 4. Write a function that calculates the molecular weight of each amino acid sequence.
## 5. Write a function that computes the GC-content of each DNA sequence.

#-- In the MAIN part of the script --#
## 6. Add two new columns to the penguin DataFrame: (1) molecular weight and (2) GC content.
## 7. Call your functions from step 3 (or step 2) and step 4 and fill in the new columns in the DataFrame.
## 8. Plot a bar-chart of adult body mass per species. In your description of the graph, provide text that answers these questions: 
#       a. What is the smallest penguin species? 
#       b. What is the geographical range of this species?
## 9. Plot a graph that shows the molecular weight as a function of GC content. 
## 10. Write the entire DataFrame to a new CSV file that includes your new columns.
## 11. BONUS: What other visualizations, functions or tasks would you do with this dataset? Add something interesting for fun. (0.5 additional points if your total score is < 15).

#-- Additional Instructions (points will be deducted if these instructions are not heeded) --#
## ** Do all of this in a Jupyter notebook and push it to a GitHub repository.
## ** Your repository should not contain any files other than those associated with this assignment. 
## ** Read all comments carefully and answer the questions by including information in your Jupyter notebook.
## ** Document all of your code (and Dr. X's code) very thoroughly so that it is clear what you did.
## ** Include in-code comments to document functions and markdown comments to document rationale.
## ** Be sure to cite (by providing URLs or other appropriate citations) information appropriately in your documented notebook.
## ** Commit and push your completed work in the Jupyter notebook to your repository.
## ** Submit the URL to your git repository via Canvas by the end of the day on May 7, 2025.

#-- Disclaimer --#
## Not all of these tasks have been covered in class and you will have to use online resources to find out how to do some of these tasks.


######################## Python Translate Script ########################

## Here's the start of our Python script. Thanks for completing it for me! - Dr. X
## IMPORTANT: install BioPython so that this will work

from Bio import SeqIO
from Bio.Data import CodonTable
import pandas as pd

#%%%%%%%%%%%%%%%#
### FUNCTIONS ###
#%%%%%%%%%%%%%%%#

## 1 ##
####### GET SEQUENCES FUNCTION ########
## Dr. X: this gets sequences 
## Please properly document this function in the Jupyter notebook 
## Your descriptions of all functions should contain information about what the function does,
## as well as information about the return types and arguments.
def get_sequences_from_file(fasta_fn):
    sequence_data_dict = {}
    for record in SeqIO.parse(fasta_fn, "fasta"):
        description = record.description.split()
        species_name = description[1] + " " + description[2]
        sequence_data_dict[species_name] = record.seq
    return(sequence_data_dict)

## 2 ##
####### YOUR STRING-TRANSLATE FUNCTION ########
## Write a function that translates sequences
## All sequences start at codon position 1
## Complete a function that translates using a loop over the string of nucleotides
## Here is  some pseudo-code and suggestions
## feel free to change the function and variable names
# def translate_function(string_nucleotides): 
#     mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"] # this should work using BioPython (be sure to check what this returns)
#     for-loop through every 3rd position in string_nucleotides to get the codon using range subsets
#         # IMPORTANT: if the sequence has a stop codon at the end, you should leave it off
#         # this is how you can retrieve the amino acid: mito_table.forward_table[codon]
#         add the aa to aa_seq_string
#     return(aa_seq_string)

## 3 ##
####### YOUR ALTERNATIVE FUNCTION ########
## Is there a better way to write the translation function? (Hint: yes there is.) 
## Perhaps using available BioPython library utilities?
## Please also write this function.


## 4 ##
####### YOUR COUNT AA ANALYSIS FUNCTION ########
## Write a function that calculates the molecular weight of each amino acid sequence.
## For this, you can use some BioPython functions. I think you can use the ProtParam module.
## For more info, check this out: http://biopython.org/wiki/ProtParam
## So you should import the following before defining your function:
from Bio.SeqUtils.ProtParam import ProteinAnalysis
# def compute_molecular_weight(aa_seq):
#     # I think the ProtParam functions may require aa_seq to be a string.
#     # It may not work if the amino acid sequence has stop codons.
#     run the ProteinAnalysis() function on aa_seq
#	  return the molecular weight

## 5 ##
####### YOUR GC CONTENT ANALYSIS FUNCTION ########
## Write a function that calculates the GC-content (proportion of "G" and "C") of each DNA sequence and returns this value.


#%%%%%%%%%%%%%%#
###   MAIN   ###
#%%%%%%%%%%%%%%#

cytb_seqs = get_sequences_from_file("penguins_cytb.fasta") 

penguins_df = pd.read_csv("penguins_mass.csv") # Includes only data for body mass 
species_list = list(penguins_df.species)

## Suggestions for steps 6 and 7
## Note there are other ways to do this that will be acceptable, as long as they are documented and work

## 6 ## 
## Add two new columns to the penguin DataFrame: (1) molecular weight and (2) GC content.
## Set the value to 'NaN' to indicate that these cells are currently empty.

## 7 ##
## Write a for-loop that translates each sequence and also gets molecular weight and computes the GC content
## of each translated sequence and adds those data to DataFrame
# for key, value in cytb_seqs.items():
#     aa_seq = translate_function(value) # whichever translation function you implement
#     get the molecular weight of aa_seq
#     get the GC content of the DNA sequence
#     fill in empty cells in DF that you created above

## 8 ##
## Plot a bar-chart of the mass with the x-axes labeled with species names.
## *Q1* What is the smallest penguin species? 
## *Q2* What is the geographical range of this species?

## 9 ##
## Plot a visualization of the molecular weight (y-axis) as a function of GC-content (x-axis).

## 10 ##
## Save the new DataFrame to a file called "penguins_mass_cytb.csv"

## 11 - BONUS ##
## What else can we do with this dataset in Python? 
## Add functions or anything that might be interesting and fun. (optional)

