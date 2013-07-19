import dssp_parse as dd

residue_max_acc = {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
                   'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0, \
                   'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0, \
                   'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0, \
                   'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0}
                   
def get_dssp_dat(name):
  dd_ob = dd.DSSPData()
  dd_ob.parseDSSP(name)
  AAs = dd_ob.getAAs()
  ACC = dd_ob.getACC()
  
  return(AAs, ACC)

def write_rsa_file(AAs, ACC, out_file):

  out = open(out_file, 'w')
  
  for site, acc in enumerate(ACC):
    if AAs[site] in residue_max_acc:
      line = str(AAs[site]) + ',' + str(float(acc)/residue_max_acc[AAs[site]]) + '\n'
      out.write(line)
      
  out.close()
  
  return 0