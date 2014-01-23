# This file reads all data for the PDB structures into R.
# Amir Shahmoradi, Wednesday  4:09 PM, Jan 23, 2014, ICMB, UT Austin
# Amir Shahmoradi, Wednesday  8:47 PM, Jan 20, 2014, ICMB, UT Austin
# Amir Shahmoradi, Wednesday 11:01 PM, Jan 15, 2014, ICMB, UT Austin
# Amir Shahmoradi, Friday 8:27 PM, Nov 8, 2013

#INPUT FILES:

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

entropy_HP = read.table('entropies/HP.entropy', header=T)
rmsf_1RD8 = read.table('molecular_dynamics/Amber/postproc/rmsf/1RD8_AB_Cpdb_CA.rmsf',header=T,comment.char="")
map_1RD8 = read.table('molecular_dynamics/map/HP_1RD8.map',header=T)
rsa_1RD8 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1RD8_sum.rsa',header=T)
dihedral_1RD8 = read.table('molecular_dynamics/Amber/postproc/dihedrals/1RD8_dihedrals.txt',header=T)
cn13_1RD8 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1RD8_sum.cn13',header=T)
wcn_1RD8 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1RD8_sum.wcn',header=T)
entropy_1RD8 = entropy_HP$entropy[!is.na(map_1RD8$pdb_pos)]
entropy_1RD8 = as.data.frame(entropy_1RD8)

data_1RD8 = data.frame(protein="1RD8", res_num=dihedral_1RD8$res_Num, res_name=dihedral_1RD8$res_name,
                       entropy=entropy_1RD8$entropy_1RD8, 
                       rsa_cr=rsa_1RD8$CRYSTAL_RSA, rsa_avg_md=rsa_1RD8$MEAN_RSA, rsa_var_md=rsa_1RD8$VAR_RSA, 
                       rmsf_avg_md=rmsf_1RD8$AvgRMSD, rmsf_std_md=rmsf_1RD8$Stdev, 
                       phi_var_md=dihedral_1RD8$var_phi, psi_var_md=dihedral_1RD8$var_psi, chi1_var_md=dihedral_1RD8$var_chi1,
                       cn13_cr=cn13_1RD8$CRYSTAL_CN, cn13_avg_md=cn13_1RD8$MEAN_CN, cn13_var_md=cn13_1RD8$VAR_CN,
                       wcn_cr=wcn_1RD8$CRYSTAL_WCN, wcn_avg_md=wcn_1RD8$MEAN_WCN, wcn_var_md=wcn_1RD8$VAR_WCN)

entropy_DPH = read.table('entropies/DPH.entropy', header=T)
rmsf_2JLY = read.table('molecular_dynamics/Amber/postproc/rmsf/2JLY_A_Cpdb_CA.rmsf',header=T,comment.char="")
map_2JLY = read.table('molecular_dynamics/map/DPH_2JLY_A.map',header=T)
rsa_2JLY = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2JLY_sum.rsa',header=T)
dihedral_2JLY = read.table('molecular_dynamics/Amber/postproc/dihedrals/2JLY_dihedrals.txt',header=T)
cn13_2JLY = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2JLY_sum.cn13',header=T)
wcn_2JLY = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2JLY_sum.wcn',header=T)
entropy_2JLY = entropy_DPH$entropy[!is.na(map_2JLY$pdb_pos)]
entropy_2JLY = as.data.frame(entropy_2JLY)

data_2JLY = data.frame(protein="2JLY", res_num=dihedral_2JLY$res_Num, res_name=dihedral_2JLY$res_name,
                       entropy=entropy_2JLY$entropy_2JLY, 
                       rsa_cr=rsa_2JLY$CRYSTAL_RSA, rsa_avg_md=rsa_2JLY$MEAN_RSA, rsa_var_md=rsa_2JLY$VAR_RSA, 
                       rmsf_avg_md=rmsf_2JLY$AvgRMSD, rmsf_std_md=rmsf_2JLY$Stdev, 
                       phi_var_md=dihedral_2JLY$var_phi, psi_var_md=dihedral_2JLY$var_psi, chi1_var_md=dihedral_2JLY$var_chi1,
                       cn13_cr=cn13_2JLY$CRYSTAL_CN, cn13_avg_md=cn13_2JLY$MEAN_CN, cn13_var_md=cn13_2JLY$VAR_CN,
                       wcn_cr=wcn_2JLY$CRYSTAL_WCN, wcn_avg_md=wcn_2JLY$MEAN_WCN, wcn_var_md=wcn_2JLY$VAR_WCN)
					   
entropy_JEHN = read.table('entropies/JEHN.entropy', header=T)
rmsf_2Z83 = read.table('molecular_dynamics/Amber/postproc/rmsf/2Z83_rmsf_ref_inpcrd.txt',header=T,comment.char="")
map_2Z83 = read.table('molecular_dynamics/map/JEHN_2Z83.map',header=T)
rsa_2Z83 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2Z83_sum.rsa',header=T)
dihedral_2Z83 = read.table('molecular_dynamics/Amber/postproc/dihedrals/2Z83_dihedrals.txt',header=T)
cn13_2Z83 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2Z83_sum.cn13',header=T)
wcn_2Z83 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2Z83_sum.wcn',header=T)
entropy_2Z83 = entropy_JEHN$entropy[!is.na(map_2Z83$pdb_pos)]
entropy_2Z83 = as.data.frame(entropy_2Z83)

data_2Z83 = data.frame(protein="2Z83", res_num=dihedral_2Z83$res_Num, res_name=dihedral_2Z83$res_name,
                       entropy=entropy_2Z83$entropy_2Z83, 
                       rsa_cr=rsa_2Z83$CRYSTAL_RSA, rsa_avg_md=rsa_2Z83$MEAN_RSA, rsa_var_md=rsa_2Z83$VAR_RSA, 
                       rmsf_avg_md=rmsf_2Z83$AvgRMSD, rmsf_std_md=rmsf_2Z83$Stdev, 
                       phi_var_md=dihedral_2Z83$var_phi, psi_var_md=dihedral_2Z83$var_psi, chi1_var_md=dihedral_2Z83$var_chi1,
                       cn13_cr=cn13_2Z83$CRYSTAL_CN, cn13_avg_md=cn13_2Z83$MEAN_CN, cn13_var_md=cn13_2Z83$VAR_CN,
                       wcn_cr=wcn_2Z83$CRYSTAL_WCN, wcn_avg_md=wcn_2Z83$MEAN_WCN, wcn_var_md=wcn_2Z83$VAR_WCN)

entropy_HCP = read.table('entropies/HCP.entropy', header=T)
rmsf_3GOL = read.table('molecular_dynamics/Amber/postproc/rmsf/3GOL_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rmsf_3GSZ = read.table('molecular_dynamics/Amber/postproc/rmsf/3GSZ_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rmsf_3I5K = read.table('molecular_dynamics/Amber/postproc/rmsf/3I5K_rmsf_ref_inpcrd.txt',header=T,comment.char="")
map_3GOL = read.table('molecular_dynamics/map/HCP_3GOL.map',header=T)
map_3GSZ = read.table('molecular_dynamics/map/HCP_3GSZ.map',header=T)
map_3I5K = read.table('molecular_dynamics/map/HCP_3I5K.map',header=T)
rsa_3GOL = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/3GOL_sum.rsa',header=T)
rsa_3GSZ = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/3GSZ_sum.rsa',header=T)
rsa_3I5K = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/3I5K_sum.rsa',header=T)
dihedral_3GOL = read.table('molecular_dynamics/Amber/postproc/dihedrals/3GOL_dihedrals.txt',header=T)
dihedral_3GSZ = read.table('molecular_dynamics/Amber/postproc/dihedrals/3GSZ_dihedrals.txt',header=T)
dihedral_3I5K = read.table('molecular_dynamics/Amber/postproc/dihedrals/3I5K_dihedrals.txt',header=T)
cn13_3GOL = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/3GOL_sum.cn13',header=T)
cn13_3GSZ = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/3GSZ_sum.cn13',header=T)
cn13_3I5K = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/3I5K_sum.cn13',header=T)
wcn_3GOL = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/3GOL_sum.wcn',header=T)
wcn_3GSZ = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/3GSZ_sum.wcn',header=T)
wcn_3I5K = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/3I5K_sum.wcn',header=T)
entropy_3GOL = entropy_HCP$entropy[!is.na(map_3GOL$pdb_pos)]
entropy_3GOL = as.data.frame(entropy_3GOL)
entropy_3GSZ = entropy_HCP$entropy[!is.na(map_3GSZ$pdb_pos)]
entropy_3GSZ = as.data.frame(entropy_3GSZ)
entropy_3I5K = entropy_HCP$entropy[!is.na(map_3I5K$pdb_pos)]
entropy_3I5K = as.data.frame(entropy_3I5K)

data_3GOL = data.frame(protein="3GOL", res_num=dihedral_3GOL$res_Num, res_name=dihedral_3GOL$res_name,
                       entropy=entropy_3GOL$entropy_3GOL, 
                       rsa_cr=rsa_3GOL$CRYSTAL_RSA, rsa_avg_md=rsa_3GOL$MEAN_RSA, rsa_var_md=rsa_3GOL$VAR_RSA, 
                       rmsf_avg_md=rmsf_3GOL$AvgRMSD, rmsf_std_md=rmsf_3GOL$Stdev, 
                       phi_var_md=dihedral_3GOL$var_phi, psi_var_md=dihedral_3GOL$var_psi, chi1_var_md=dihedral_3GOL$var_chi1,
                       cn13_cr=cn13_3GOL$CRYSTAL_CN, cn13_avg_md=cn13_3GOL$MEAN_CN, cn13_var_md=cn13_3GOL$VAR_CN,
                       wcn_cr=wcn_3GOL$CRYSTAL_WCN, wcn_avg_md=wcn_3GOL$MEAN_WCN, wcn_var_md=wcn_3GOL$VAR_WCN)

data_3GSZ = data.frame(protein="3GSZ", res_num=dihedral_3GSZ$res_Num, res_name=dihedral_3GSZ$res_name,
                       entropy=entropy_3GSZ$entropy_3GSZ, 
                       rsa_cr=rsa_3GSZ$CRYSTAL_RSA, rsa_avg_md=rsa_3GSZ$MEAN_RSA, rsa_var_md=rsa_3GSZ$VAR_RSA, 
                       rmsf_avg_md=rmsf_3GSZ$AvgRMSD, rmsf_std_md=rmsf_3GSZ$Stdev, 
                       phi_var_md=dihedral_3GSZ$var_phi, psi_var_md=dihedral_3GSZ$var_psi, chi1_var_md=dihedral_3GSZ$var_chi1,
                       cn13_cr=cn13_3GSZ$CRYSTAL_CN, cn13_avg_md=cn13_3GSZ$MEAN_CN, cn13_var_md=cn13_3GSZ$VAR_CN,
                       wcn_cr=wcn_3GSZ$CRYSTAL_WCN, wcn_avg_md=wcn_3GSZ$MEAN_WCN, wcn_var_md=wcn_3GSZ$VAR_WCN)

data_3I5K = data.frame(protein="3I5K", res_num=dihedral_3I5K$res_Num, res_name=dihedral_3I5K$res_name,
                       entropy=entropy_3I5K$entropy_3I5K, 
                       rsa_cr=rsa_3I5K$CRYSTAL_RSA, rsa_avg_md=rsa_3I5K$MEAN_RSA, rsa_var_md=rsa_3I5K$VAR_RSA, 
                       rmsf_avg_md=rmsf_3I5K$AvgRMSD, rmsf_std_md=rmsf_3I5K$Stdev, 
                       phi_var_md=dihedral_3I5K$var_phi, psi_var_md=dihedral_3I5K$var_psi, chi1_var_md=dihedral_3I5K$var_chi1,
                       cn13_cr=cn13_3I5K$CRYSTAL_CN, cn13_avg_md=cn13_3I5K$MEAN_CN, cn13_var_md=cn13_3I5K$VAR_CN,
                       wcn_cr=wcn_3I5K$CRYSTAL_WCN, wcn_avg_md=wcn_3I5K$MEAN_WCN, wcn_var_md=wcn_3I5K$VAR_WCN)

entropy_RVFVNP = read.table('entropies/RVFVNP.entropy', header=T)
rmsf_3LYF = read.table('molecular_dynamics/Amber/postproc/rmsf/3LYF_rmsf_ref_inpcrd.txt',header=T,comment.char="")
map_3LYF = read.table('molecular_dynamics/map/RVFVNP_3LYF.map',header=T)
rsa_3LYF = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/3LYF_sum.rsa',header=T)
dihedral_3LYF = read.table('molecular_dynamics/Amber/postproc/dihedrals/3LYF_dihedrals.txt',header=T)
cn13_3LYF = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/3LYF_sum.cn13',header=T)
wcn_3LYF = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/3LYF_sum.wcn',header=T)
entropy_3LYF = entropy_RVFVNP$entropy[!is.na(map_3LYF$pdb_pos)]
entropy_3LYF = as.data.frame(entropy_3LYF)

data_3LYF = data.frame(protein="3LYF", res_num=dihedral_3LYF$res_Num, res_name=dihedral_3LYF$res_name,
                       entropy=entropy_3LYF$entropy_3LYF,
                       rsa_cr=rsa_3LYF$CRYSTAL_RSA, rsa_avg_md=rsa_3LYF$MEAN_RSA, rsa_var_md=rsa_3LYF$VAR_RSA, 
                       rmsf_avg_md=rmsf_3LYF$AvgRMSD, rmsf_std_md=rmsf_3LYF$Stdev, 
                       phi_var_md=dihedral_3LYF$var_phi, psi_var_md=dihedral_3LYF$var_psi, chi1_var_md=dihedral_3LYF$var_chi1,
                       cn13_cr=cn13_3LYF$CRYSTAL_CN, cn13_avg_md=cn13_3LYF$MEAN_CN, cn13_var_md=cn13_3LYF$VAR_CN,
                       wcn_cr=wcn_3LYF$CRYSTAL_WCN, wcn_avg_md=wcn_3LYF$MEAN_WCN, wcn_var_md=wcn_3LYF$VAR_WCN)
					   
entropy_CCHFN = read.table('entropies/CCHFN.entropy', header=T)
rmsf_4AQF_B = read.table('molecular_dynamics/Amber/postproc/rmsf/4AQF_B_rmsf_ref_inpcrd.txt',header=T,comment.char="")
map_4AQF_B = read.table('molecular_dynamics/map/CCHFN_4AQF_B.map',header=T)
rsa_4AQF_B = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/4AQF_B_sum.rsa',header=T)
dihedral_4AQF_B = read.table('molecular_dynamics/Amber/postproc/dihedrals/4AQF_B_dihedrals.txt',header=T)
cn13_4AQF_B = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/4AQF_B_sum.cn13',header=T)
wcn_4AQF_B = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/4AQF_B_sum.wcn',header=T)
entropy_4AQF_B = entropy_CCHFN$entropy[!is.na(map_4AQF_B$pdb_pos)]
entropy_4AQF_B = as.data.frame(entropy_4AQF_B)

data_4AQF_B = data.frame(protein="4AQF_B", res_num=dihedral_4AQF_B$res_Num, res_name=dihedral_4AQF_B$res_name,
                       entropy=entropy_4AQF_B$entropy_4AQF_B, 
                       rsa_cr=rsa_4AQF_B$CRYSTAL_RSA, rsa_avg_md=rsa_4AQF_B$MEAN_RSA, rsa_var_md=rsa_4AQF_B$VAR_RSA, 
                       rmsf_avg_md=rmsf_4AQF_B$AvgRMSD, rmsf_std_md=rmsf_4AQF_B$Stdev, 
                       phi_var_md=dihedral_4AQF_B$var_phi, psi_var_md=dihedral_4AQF_B$var_psi, chi1_var_md=dihedral_4AQF_B$var_chi1,
                       cn13_cr=cn13_4AQF_B$CRYSTAL_CN, cn13_avg_md=cn13_4AQF_B$MEAN_CN, cn13_var_md=cn13_4AQF_B$VAR_CN,
                       wcn_cr=wcn_4AQF_B$CRYSTAL_WCN, wcn_avg_md=wcn_4AQF_B$MEAN_WCN, wcn_var_md=wcn_4AQF_B$VAR_WCN)
					   
entropy_MRNABD = read.table('entropies/MRNABD.entropy', header=T)
rmsf_4GHA_A = read.table('molecular_dynamics/Amber/postproc/rmsf/4GHA_A_rmsf_ref_inpcrd.txt',header=T,comment.char="")
map_4GHA_A = read.table('molecular_dynamics/map/MRNABD_4GHA_A.map',header=T)
rsa_4GHA_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/4GHA_A_sum.rsa',header=T)
dihedral_4GHA_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/4GHA_A_dihedrals.txt',header=T)
cn13_4GHA_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/4GHA_A_sum.cn13',header=T)
wcn_4GHA_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/4GHA_A_sum.wcn',header=T)
entropy_4GHA_A = entropy_MRNABD$entropy[!is.na(map_4GHA_A$pdb_pos)]
entropy_4GHA_A = as.data.frame(entropy_4GHA_A)

data_4GHA_A = data.frame(protein="4GHA_A", res_num=dihedral_4GHA_A$res_Num, res_name=dihedral_4GHA_A$res_name,
                       entropy=entropy_4GHA_A$entropy_4GHA_A, 
                       rsa_cr=rsa_4GHA_A$CRYSTAL_RSA, rsa_avg_md=rsa_4GHA_A$MEAN_RSA, rsa_var_md=rsa_4GHA_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_4GHA_A$AvgRMSD, rmsf_std_md=rmsf_4GHA_A$Stdev, 
                       phi_var_md=dihedral_4GHA_A$var_phi, psi_var_md=dihedral_4GHA_A$var_psi, chi1_var_md=dihedral_4GHA_A$var_chi1,
                       cn13_cr=cn13_4GHA_A$CRYSTAL_CN, cn13_avg_md=cn13_4GHA_A$MEAN_CN, cn13_var_md=cn13_4GHA_A$VAR_CN,
                       wcn_cr=wcn_4GHA_A$CRYSTAL_WCN, wcn_avg_md=wcn_4GHA_A$MEAN_WCN, wcn_var_md=wcn_4GHA_A$VAR_WCN)

entropy_INP = read.table('entropies/INP.entropy', header=T)
rmsf_4IRY = read.table('molecular_dynamics/Amber/postproc/rmsf/4IRY_rmsf_ref_inpcrd.txt',header=T,comment.char="")
map_4IRY = read.table('molecular_dynamics/map/INP_4IRY.map',header=T)
rsa_4IRY = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/4IRY_sum.rsa',header=T)
dihedral_4IRY = read.table('molecular_dynamics/Amber/postproc/dihedrals/4IRY_dihedrals.txt',header=T)
cn13_4IRY = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/4IRY_sum.cn13',header=T)
wcn_4IRY = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/4IRY_sum.wcn',header=T)
entropy_4IRY = entropy_INP$entropy[!is.na(map_4IRY$pdb_pos)]
entropy_4IRY = as.data.frame(entropy_4IRY)

data_4IRY = data.frame(protein="4IRY", res_num=dihedral_4IRY$res_Num, res_name=dihedral_4IRY$res_name,
                       entropy=entropy_4IRY$entropy_4IRY, 
                       rsa_cr=rsa_4IRY$CRYSTAL_RSA, rsa_avg_md=rsa_4IRY$MEAN_RSA, rsa_var_md=rsa_4IRY$VAR_RSA, 
                       rmsf_avg_md=rmsf_4IRY$AvgRMSD, rmsf_std_md=rmsf_4IRY$Stdev, 
                       phi_var_md=dihedral_4IRY$var_phi, psi_var_md=dihedral_4IRY$var_psi, chi1_var_md=dihedral_4IRY$var_chi1,
                       cn13_cr=cn13_4IRY$CRYSTAL_CN, cn13_avg_md=cn13_4IRY$MEAN_CN, cn13_var_md=cn13_4IRY$VAR_CN,
                       wcn_cr=wcn_4IRY$CRYSTAL_WCN, wcn_avg_md=wcn_4IRY$MEAN_WCN, wcn_var_md=wcn_4IRY$VAR_WCN)

entropy_WNPB = read.table('entropies/WNPB.entropy', header=T)
rmsf_2FP7_B = read.table('molecular_dynamics/Amber/postproc/rmsf/2FP7_B_rmsf_CA_ref_pdb.txt',header=T,comment.char="")
map_2FP7_B = read.table('molecular_dynamics/map/WNPB_2FP7_B.map',header=T)
rsa_2FP7_B = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2FP7_B_sum.rsa',header=T)
dihedral_2FP7_B = read.table('molecular_dynamics/Amber/postproc/dihedrals/2FP7_B_dihedrals.txt',header=T)
cn13_2FP7_B = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2FP7_B_sum.cn13',header=T)
wcn_2FP7_B = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2FP7_B_sum.wcn',header=T)
entropy_2FP7_B = entropy_WNPB$entropy[!is.na(map_2FP7_B$pdb_pos)]
entropy_2FP7_B = as.data.frame(entropy_2FP7_B)

data_2FP7_B = data.frame(protein="2FP7_B", res_num=dihedral_2FP7_B$res_Num, res_name=dihedral_2FP7_B$res_name,
                       entropy=entropy_2FP7_B$entropy_2FP7_B, 
                       rsa_cr=rsa_2FP7_B$CRYSTAL_RSA, rsa_avg_md=rsa_2FP7_B$MEAN_RSA, rsa_var_md=rsa_2FP7_B$VAR_RSA, 
                       rmsf_avg_md=rmsf_2FP7_B$AvgRMSD, rmsf_std_md=rmsf_2FP7_B$Stdev, 
                       phi_var_md=dihedral_2FP7_B$var_phi, psi_var_md=dihedral_2FP7_B$var_psi, chi1_var_md=dihedral_2FP7_B$var_chi1,
                       cn13_cr=cn13_2FP7_B$CRYSTAL_CN, cn13_avg_md=cn13_2FP7_B$MEAN_CN, cn13_var_md=cn13_2FP7_B$VAR_CN,
                       wcn_cr=wcn_2FP7_B$CRYSTAL_WCN, wcn_avg_md=wcn_2FP7_B$MEAN_WCN, wcn_var_md=wcn_2FP7_B$VAR_WCN)

# THERMOPHILIC STRUCTURES

rmsf_1AJ8_A = read.table('molecular_dynamics/Amber/postproc/rmsf/1AJ8_A_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rsa_1AJ8_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1AJ8_A_sum.rsa',header=T)
dihedral_1AJ8_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/1AJ8_A_dihedrals.txt',header=T)
cn13_1AJ8_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1AJ8_A_sum.cn13',header=T)
wcn_1AJ8_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1AJ8_A_sum.wcn',header=T)

data_1AJ8_A = data.frame(protein="1AJ8_A", res_num=dihedral_1AJ8_A$res_Num, res_name=dihedral_1AJ8_A$res_name,
                       rsa_cr=rsa_1AJ8_A$CRYSTAL_RSA, rsa_avg_md=rsa_1AJ8_A$MEAN_RSA, rsa_var_md=rsa_1AJ8_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_1AJ8_A$AvgRMSD, rmsf_std_md=rmsf_1AJ8_A$Stdev, 
                       phi_var_md=dihedral_1AJ8_A$var_phi, psi_var_md=dihedral_1AJ8_A$var_psi, chi1_var_md=dihedral_1AJ8_A$var_chi1,
                       cn13_cr=cn13_1AJ8_A$CRYSTAL_CN, cn13_avg_md=cn13_1AJ8_A$MEAN_CN, cn13_var_md=cn13_1AJ8_A$VAR_CN,
                       wcn_cr=wcn_1AJ8_A$CRYSTAL_WCN, wcn_avg_md=wcn_1AJ8_A$MEAN_WCN, wcn_var_md=wcn_1AJ8_A$VAR_WCN)

rmsf_1AOR_A = read.table('molecular_dynamics/Amber/postproc/rmsf/1AOR_A_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rsa_1AOR_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1AOR_A_sum.rsa',header=T)
dihedral_1AOR_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/1AOR_A_dihedrals.txt',header=T)
cn13_1AOR_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1AOR_A_sum.cn13',header=T)
wcn_1AOR_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1AOR_A_sum.wcn',header=T)

data_1AOR_A = data.frame(protein="1AOR_A", res_num=dihedral_1AOR_A$res_Num, res_name=dihedral_1AOR_A$res_name,
                       rsa_cr=rsa_1AOR_A$CRYSTAL_RSA, rsa_avg_md=rsa_1AOR_A$MEAN_RSA, rsa_var_md=rsa_1AOR_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_1AOR_A$AvgRMSD, rmsf_std_md=rmsf_1AOR_A$Stdev, 
                       phi_var_md=dihedral_1AOR_A$var_phi, psi_var_md=dihedral_1AOR_A$var_psi, chi1_var_md=dihedral_1AOR_A$var_chi1,
                       cn13_cr=cn13_1AOR_A$CRYSTAL_CN, cn13_avg_md=cn13_1AOR_A$MEAN_CN, cn13_var_md=cn13_1AOR_A$VAR_CN,
                       wcn_cr=wcn_1AOR_A$CRYSTAL_WCN, wcn_avg_md=wcn_1AOR_A$MEAN_WCN, wcn_var_md=wcn_1AOR_A$VAR_WCN)

rmsf_1CTS_A = read.table('molecular_dynamics/Amber/postproc/rmsf/1CTS_A_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rsa_1CTS_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1CTS_A_sum.rsa',header=T)
dihedral_1CTS_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/1CTS_A_dihedrals.txt',header=T)
cn13_1CTS_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1CTS_A_sum.cn13',header=T)
wcn_1CTS_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1CTS_A_sum.wcn',header=T)

data_1CTS_A = data.frame(protein="1CTS_A", res_num=dihedral_1CTS_A$res_Num, res_name=dihedral_1CTS_A$res_name,
                       rsa_cr=rsa_1CTS_A$CRYSTAL_RSA, rsa_avg_md=rsa_1CTS_A$MEAN_RSA, rsa_var_md=rsa_1CTS_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_1CTS_A$AvgRMSD, rmsf_std_md=rmsf_1CTS_A$Stdev, 
                       phi_var_md=dihedral_1CTS_A$var_phi, psi_var_md=dihedral_1CTS_A$var_psi, chi1_var_md=dihedral_1CTS_A$var_chi1,
                       cn13_cr=cn13_1CTS_A$CRYSTAL_CN, cn13_avg_md=cn13_1CTS_A$MEAN_CN, cn13_var_md=cn13_1CTS_A$VAR_CN,
                       wcn_cr=wcn_1CTS_A$CRYSTAL_WCN, wcn_avg_md=wcn_1CTS_A$MEAN_WCN, wcn_var_md=wcn_1CTS_A$VAR_WCN)

rmsf_1MP9_A = read.table('molecular_dynamics/Amber/postproc/rmsf/1MP9_A_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rsa_1MP9_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1MP9_A_sum.rsa',header=T)
dihedral_1MP9_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/1MP9_A_dihedrals.txt',header=T)
cn13_1MP9_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1MP9_A_sum.cn13',header=T)
wcn_1MP9_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1MP9_A_sum.wcn',header=T)

data_1MP9_A = data.frame(protein="1MP9_A", res_num=dihedral_1MP9_A$res_Num, res_name=dihedral_1MP9_A$res_name,
                       rsa_cr=rsa_1MP9_A$CRYSTAL_RSA, rsa_avg_md=rsa_1MP9_A$MEAN_RSA, rsa_var_md=rsa_1MP9_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_1MP9_A$AvgRMSD, rmsf_std_md=rmsf_1MP9_A$Stdev, 
                       phi_var_md=dihedral_1MP9_A$var_phi, psi_var_md=dihedral_1MP9_A$var_psi, chi1_var_md=dihedral_1MP9_A$var_chi1,
                       cn13_cr=cn13_1MP9_A$CRYSTAL_CN, cn13_avg_md=cn13_1MP9_A$MEAN_CN, cn13_var_md=cn13_1MP9_A$VAR_CN,
                       wcn_cr=wcn_1MP9_A$CRYSTAL_WCN, wcn_avg_md=wcn_1MP9_A$MEAN_WCN, wcn_var_md=wcn_1MP9_A$VAR_WCN)

# DIFFERENT TEMPERATURE

rmsf_2JLY_temp_50 = read.table('molecular_dynamics/Amber/postproc/rmsf/2JLY_temp_50_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rsa_2JLY_temp_50 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2JLY_temp_50_sum.rsa',header=T)
dihedral_2JLY_temp_50 = read.table('molecular_dynamics/Amber/postproc/dihedrals/2JLY_temp_50_dihedrals.txt',header=T)
cn13_2JLY_temp_50 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2JLY_temp_50_sum.cn13',header=T)
wcn_2JLY_temp_50 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2JLY_temp_50_sum.wcn',header=T)

data_2JLY_temp_50 = data.frame(protein="2JLY_temp_50", res_num=dihedral_2JLY$res_Num, res_name=dihedral_2JLY$res_name,
                       entropy=entropy_2JLY$entropy_2JLY, 
                       rsa_cr=rsa_2JLY_temp_50$CRYSTAL_RSA, rsa_avg_md=rsa_2JLY_temp_50$MEAN_RSA, rsa_var_md=rsa_2JLY_temp_50$VAR_RSA, 
                       rmsf_avg_md=rmsf_2JLY_temp_50$AvgRMSD, rmsf_std_md=rmsf_2JLY_temp_50$Stdev, 
                       phi_var_md=dihedral_2JLY_temp_50$var_phi, psi_var_md=dihedral_2JLY_temp_50$var_psi, chi1_var_md=dihedral_2JLY_temp_50$var_chi1,
                       cn13_cr=cn13_2JLY_temp_50$CRYSTAL_CN, cn13_avg_md=cn13_2JLY_temp_50$MEAN_CN, cn13_var_md=cn13_2JLY_temp_50$VAR_CN,
                       wcn_cr=wcn_2JLY_temp_50$CRYSTAL_WCN, wcn_avg_md=wcn_2JLY_temp_50$MEAN_WCN, wcn_var_md=wcn_2JLY_temp_50$VAR_WCN)

rmsf_2JLY_temp_100 = read.table('molecular_dynamics/Amber/postproc/rmsf/2JLY_temp_100_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rsa_2JLY_temp_100 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2JLY_temp_100_sum.rsa',header=T)
dihedral_2JLY_temp_100 = read.table('molecular_dynamics/Amber/postproc/dihedrals/2JLY_temp_100_dihedrals.txt',header=T)
cn13_2JLY_temp_100 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2JLY_temp_100_sum.cn13',header=T)
wcn_2JLY_temp_100 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2JLY_temp_100_sum.wcn',header=T)

data_2JLY_temp_100 = data.frame(protein="2JLY_temp_100", res_num=dihedral_2JLY$res_Num, res_name=dihedral_2JLY$res_name,
                       entropy=entropy_2JLY$entropy_2JLY, 
                       rsa_cr=rsa_2JLY_temp_100$CRYSTAL_RSA, rsa_avg_md=rsa_2JLY_temp_100$MEAN_RSA, rsa_var_md=rsa_2JLY_temp_100$VAR_RSA, 
                       rmsf_avg_md=rmsf_2JLY_temp_100$AvgRMSD, rmsf_std_md=rmsf_2JLY_temp_100$Stdev, 
                       phi_var_md=dihedral_2JLY_temp_100$var_phi, psi_var_md=dihedral_2JLY_temp_100$var_psi, chi1_var_md=dihedral_2JLY_temp_100$var_chi1,
                       cn13_cr=cn13_2JLY_temp_100$CRYSTAL_CN, cn13_avg_md=cn13_2JLY_temp_100$MEAN_CN, cn13_var_md=cn13_2JLY_temp_100$VAR_CN,
                       wcn_cr=wcn_2JLY_temp_100$CRYSTAL_WCN, wcn_avg_md=wcn_2JLY_temp_100$MEAN_WCN, wcn_var_md=wcn_2JLY_temp_100$VAR_WCN)

rmsf_2JLY_temp_200 = read.table('molecular_dynamics/Amber/postproc/rmsf/2JLY_temp_200_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rsa_2JLY_temp_200 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2JLY_temp_200_sum.rsa',header=T)
dihedral_2JLY_temp_200 = read.table('molecular_dynamics/Amber/postproc/dihedrals/2JLY_temp_200_dihedrals.txt',header=T)
cn13_2JLY_temp_200 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2JLY_temp_200_sum.cn13',header=T)
wcn_2JLY_temp_200 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2JLY_temp_200_sum.wcn',header=T)

data_2JLY_temp_200 = data.frame(protein="2JLY_temp_200", res_num=dihedral_2JLY$res_Num, res_name=dihedral_2JLY$res_name,
                       entropy=entropy_2JLY$entropy_2JLY, 
                       rsa_cr=rsa_2JLY_temp_200$CRYSTAL_RSA, rsa_avg_md=rsa_2JLY_temp_200$MEAN_RSA, rsa_var_md=rsa_2JLY_temp_200$VAR_RSA, 
                       rmsf_avg_md=rmsf_2JLY_temp_200$AvgRMSD, rmsf_std_md=rmsf_2JLY_temp_200$Stdev, 
                       phi_var_md=dihedral_2JLY_temp_200$var_phi, psi_var_md=dihedral_2JLY_temp_200$var_psi, chi1_var_md=dihedral_2JLY_temp_200$var_chi1,
                       cn13_cr=cn13_2JLY_temp_200$CRYSTAL_CN, cn13_avg_md=cn13_2JLY_temp_200$MEAN_CN, cn13_var_md=cn13_2JLY_temp_200$VAR_CN,
                       wcn_cr=wcn_2JLY_temp_200$CRYSTAL_WCN, wcn_avg_md=wcn_2JLY_temp_200$MEAN_WCN, wcn_var_md=wcn_2JLY_temp_200$VAR_WCN)

rmsf_2JLY_temp_450 = read.table('molecular_dynamics/Amber/postproc/rmsf/2JLY_temp_450_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rsa_2JLY_temp_450 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2JLY_temp_450_sum.rsa',header=T)
dihedral_2JLY_temp_450 = read.table('molecular_dynamics/Amber/postproc/dihedrals/2JLY_temp_450_dihedrals.txt',header=T)
cn13_2JLY_temp_450 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2JLY_temp_450_sum.cn13',header=T)
wcn_2JLY_temp_450 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2JLY_temp_450_sum.wcn',header=T)

data_2JLY_temp_450 = data.frame(protein="2JLY_temp_450", res_num=dihedral_2JLY$res_Num, res_name=dihedral_2JLY$res_name,
                       entropy=entropy_2JLY$entropy_2JLY, 
                       rsa_cr=rsa_2JLY_temp_450$CRYSTAL_RSA, rsa_avg_md=rsa_2JLY_temp_450$MEAN_RSA, rsa_var_md=rsa_2JLY_temp_450$VAR_RSA, 
                       rmsf_avg_md=rmsf_2JLY_temp_450$AvgRMSD, rmsf_std_md=rmsf_2JLY_temp_450$Stdev, 
                       phi_var_md=dihedral_2JLY_temp_450$var_phi, psi_var_md=dihedral_2JLY_temp_450$var_psi, chi1_var_md=dihedral_2JLY_temp_450$var_chi1,
                       cn13_cr=cn13_2JLY_temp_450$CRYSTAL_CN, cn13_avg_md=cn13_2JLY_temp_450$MEAN_CN, cn13_var_md=cn13_2JLY_temp_450$VAR_CN,
                       wcn_cr=wcn_2JLY_temp_450$CRYSTAL_WCN, wcn_avg_md=wcn_2JLY_temp_450$MEAN_WCN, wcn_var_md=wcn_2JLY_temp_450$VAR_WCN)

rmsf_1AOR_A_temp_373 = read.table('molecular_dynamics/Amber/postproc/rmsf/1AOR_A_temp_373_rmsf_ref_inpcrd.txt',header=T,comment.char="")
rsa_1AOR_A_temp_373 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1AOR_A_temp_373_sum.rsa',header=T)
dihedral_1AOR_A_temp_373 = read.table('molecular_dynamics/Amber/postproc/dihedrals/1AOR_A_temp_373_dihedrals.txt',header=T)
cn13_1AOR_A_temp_373 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1AOR_A_temp_373_sum.cn13',header=T)
wcn_1AOR_A_temp_373 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1AOR_A_temp_373_sum.wcn',header=T)

data_1AOR_A_temp_373 = data.frame(protein="1AOR_A", res_num=dihedral_1AOR_A$res_Num, res_name=dihedral_1AOR_A$res_name,
                       rsa_cr=rsa_1AOR_A_temp_373$CRYSTAL_RSA, rsa_avg_md=rsa_1AOR_A_temp_373$MEAN_RSA, rsa_var_md=rsa_1AOR_A_temp_373$VAR_RSA, 
                       rmsf_avg_md=rmsf_1AOR_A_temp_373$AvgRMSD, rmsf_std_md=rmsf_1AOR_A_temp_373$Stdev, 
                       phi_var_md=dihedral_1AOR_A_temp_373$var_phi, psi_var_md=dihedral_1AOR_A_temp_373$var_psi, chi1_var_md=dihedral_1AOR_A_temp_373$var_chi1,
                       cn13_cr=cn13_1AOR_A_temp_373$CRYSTAL_CN, cn13_avg_md=cn13_1AOR_A_temp_373$MEAN_CN, cn13_var_md=cn13_1AOR_A_temp_373$VAR_CN,
                       wcn_cr=wcn_1AOR_A_temp_373$CRYSTAL_WCN, wcn_avg_md=wcn_1AOR_A_temp_373$MEAN_WCN, wcn_var_md=wcn_1AOR_A_temp_373$VAR_WCN)

