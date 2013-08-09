#BEFORE RUNNING THIS SCRIPT!!!! in all files, do these 3 search and replace w/ nothing:      {      }     ^\s+

## 8/9/13 by SJS for the ASAP project.

import re, os, fnmatch, csv

#################################################################
def getSiteRates(post_matrix,  omegas):
	'''Create a list of rates at each position. Compute each as a weighted average of sum(postprob * w).'''

	allrates=[]
	for n in range(len(post_matrix[0])):
		rate = 0
		for w in range(0,len(omegas)):
			temp = float(postmat[w][n]) * float(omegas[w])	
			rate = rate + temp
		allrates.append(rate)
	
	return allrates
#################################################################



## 'hivpol',
names=['riftvalley',  'H1N1_HA', 'dengue_ns3', 'H1N1_NP', 'HCV_NS5B', 'marburg.vp35', 'JEV', 'westNile_ns2b', 'westnile_chainb', 'crimean_congo.Nucleoprotein']

hyout_dir = '/Users/sjspielman/Dropbox/ViralProject_SJS/HyPhy/hyout/'
numcat=5 # We fixed there to be 5 rate categories for every gene.

for name in names:
	print name
	
	file = open(hyout_dir+name+'.out')
	filestring = file.read()
	file.close()
	file = open(hyout_dir+name+'.out')
	lines=file.readlines()
	numlines=len(lines)
	file.close()

	## Retrieve the rate categories (w_1, w_2, w_3, w_4, w_5)
	omegas=[]
	for w in range(numcat+1):
		find=re.search('w_'+str(w)+'=(\d\.*\d*)', filestring)
		find_scinot = re.search('w_'+str(w)+'=\d\.*\d*e-\d+', filestring)
		if find != None:
			if find_scinot:
				omegas.append(float(0))  ## If small enough to resort to scientific notation, can round it to 0.
			else:
				omegas.append(float(find.group(1)))
				
	## Get the posterior probability matrix of postprobs for site rates. This is the THIRD matrix in the file (as in, the last one at the bottom)
	postmat=[]
	for n in range(numlines-numcat, numlines):
		entry=''
		probs=[]
		probs_raw=lines[n]
		for x in range(0, len(probs_raw)):
			if probs_raw[x]!=',' and probs_raw[x]!='\n':
				entry=entry+probs_raw[x]
			else:
				probs.append(float(entry))
				entry=''
		postmat.append(probs)

	## Compute rates at each site by weighting omegas with post probs
	rates = getSiteRates(postmat, omegas)
	
	## Save to file
	outfile = '/Users/sjspielman/Dropbox/ViralProject_SJS/siterates_REL/'+name+'.txt'
	outhandle = open(outfile, 'w')
	sitecount=1
	for r in rates:
		outhandle.write(str(sitecount)+'\t'+str(r)+'\n')
		sitecount+=1
	outhandle.close()
	
	