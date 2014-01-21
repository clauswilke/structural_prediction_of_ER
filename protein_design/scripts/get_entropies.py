#!usr/local/bin/python
import math, string, re
import numpy as np

#	Written/Last Updated by Eleisha Jackson on August 8, 2013
#	Description: This script calculates the entropies a list of PDBS. This script must be in the same
#	directory as the files with the alignments needed to calculate entropy.

#This is a dictionary that maps the three-letter amino acid to to its one letter code
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',                   
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',                   
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',                   
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }          

PDBS = ["1RD8", "2JLY", "2Z83", "3GOL", "3GSZ", "3I5K", "3LYF", "4AQF_B", "4GHA_A", "4IRY"] #List of PDBS hard-coded that you want to calculate entropies for
chain = 'X' #For these proteins they each have the same chain, this can changed if that is not the case. 
#chain = 'B'
#PDBS = ["2FP7"]
temp = 0.6 #Temperature for the designed proteins that are being analyzed.

#This function takes the amino acid count data for the designed and corresponding natural sequences and returns a list of lists with all the AA data at each site. 
def get_AA_distribution(AA_counts): 
    protein_distribution = []
    transformed_protein_distribution = []
    num_AA = 0
    protein_distribution = AA_counts #Set the protein distribution as the counts used as niput
    for site in protein_distribution: #For each site in the list with the lists of amino acid counts
        new_site = site
        num_AA = sum(new_site)
        aa_probs = []
        for count in new_site: #This looks turns the counts in the site into frequencies 
            if count == 0:
                prob = float(1)/float(num_AA + 20)
            else:
                prob = float(count + 1)/float(num_AA + 20)
            aa_probs.append(prob)
        transformed_protein_distribution.append(aa_probs)  #Appends the list of frequencies for each site to a giant list
    #m,n = np.array(transformed_protein_distribution).shape 
    return transformed_protein_distribution	 

#This takes a file with site amino acid count data and then calculates the entropy for sites. 
def calculate_entropy(AA_counts):
    probs = get_AA_distribution(AA_counts) #Get the amino acid site data
    entropy_values = []
    entropy_number = 0
    probs_array = np.array(probs)
    num_residues,num_AA = probs_array.shape
        
    for i in xrange(0, num_residues): #For each site
        probs_values = probs_array[i]
        prob_sum = sum(probs_values)
        entropy_number = 0
        for j in xrange(0,20): #Calculate the entropy (Can look up this formula, just the native entropy)
            value = (float(probs_values[j])*np.log(float(probs_values[j])))
            entropy_number = entropy_number + value
        entropy_values.append(-entropy_number) #Append entropy value to entropy at sites. 
    return entropy_values

#This is a helper function that takes a string and formats it with commas between each character.
#Useful for splitting strings into lists (can then use re to split by comma)
def dump_csv_line(line):
    new_line = ""
    size = len(line)
    for x in range(0,size):
        new_line = new_line + str(line[x])
        if(x != size-1):
            new_line = new_line + ","

    new_line += "\n"
    return new_line

#This functions takes a file with sequences and then returns a lists of lists. Each list in the list of lists is just
#a list of AAs in each sequence of the file. 
def get_AA_lists(sequence_filename):	
	file = open(sequence_filename, "r")
	seq_lists = []
	new_sequences = []
	sequences = file.readlines()
	for seq in sequences:
		new_seq = dump_csv_line(seq) #Split each sequence in the alignment by commas between the amino acids
		new_sequences.append(new_seq)
	for seq in new_sequences:
		AA_list = re.split(",",seq) #Make an array of 
		AA_list.pop()
		seq_lists.append(AA_list)
	return seq_lists
	
#This file takes as its input a list of lists where each list is made of amino acids. It then calculates the number of times 
#that an amino acid is seen as a site. Each element each list is a site in a sequence. 
#Note: Each list in the input of list of lists should have the same length. 
def get_AA_counts(seq_lists):	
	counter = 0
	aaSum = 0			  
	AA = resdict.values() #Gets a list of all the one letter-codes for the amino acids
	data = np.array(seq_lists)
	seq_length = len(seq_lists[0])
	aaCount_lists = []
	n,m = data.shape 
	for i in range(0, m): #Loops through for each column in the array (Loops over all sites in the alignment)
		aaCount = [] #Makes a list to store the frequency data for each site
		try:
			aaList = list(data[1:,i]) #Gets the amino acids counts for that site
		except IndexError:
			print "--------------"
			print "Something wrong with list of AAs"
			print i
			print data.ndim
			print data.shape
			print "--------------"
			quit()	
		for aa in AA: 
			try:
				aaCount.append(aaList.count(aa)) #Counts the frequencies for each amino acid at this site
			except ValueError,IndexError:
				print "-----------------------"
				print "Weird AA Key!"
				print aaList
				print "-----------------------"
				quit()
		aaCount_lists.append(aaCount) #Appends the list of frequency data for this site to a list 
	return aaCount_lists #Returns a list of lists of the frequency data for each site. Each list should have length 20.

for PDB in PDBS: #For each protein in the PDB list
	file = "sequences_" + PDB + "_" + chain + "_" + str(temp) + ".csv" #Get the file where the aligned sequences are store
	print "Processing: ", file 
	data_lists = get_AA_lists(file) #Get a list of lists with the amino acids of each sequence in the alignment
	protein_site_counts = get_AA_counts(data_lists) #Get the frequency data for each site of 
	#a,b = np.array(protein_site_counts).shape
	entropies = calculate_entropy(protein_site_counts) #Calculate the entropy for site in the alignment
	for j in entropies: #This loop checks to make sure the entropies are reasonable. Should not be greater than log(n),
		if j > np.log(20): #In this case n = 20 for 20 amino acids
			print "The entropy is too high! Mathematically impossible!"
			print j
	output = PDB + "_" + chain + "_" + str(temp) + "_entropy.txt" #Creates output file to print entropies for each site
	out_file = open(output, "w")	
	out_file.write("res_num\tentropy\n")
	for i in range(0, len(entropies)): #This loop prints the entropy data for each site to the output file
		if (i == (len(entropies)-1)):
			entropy_string = str(i+1) + "\t" + str(entropies[i])
		else:
			entropy_string = str(i+1) + "\t" + str(entropies[i]) + "\n"
		out_file.write(entropy_string)
	out_file.close() #Close the output that was opened

	
	