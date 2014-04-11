# This file reads all data for the PDB structures into R.
# Amir Shahmoradi,  Saturday 10:38 PM, Feb  1, 2014, ICMB, UT Austin
# Amir Shahmoradi,  Thursday  4:09 PM, Jan 23, 2014, ICMB, UT Austin
# Amir Shahmoradi,    Monday  8:47 PM, Jan 20, 2014, ICMB, UT Austin
# Amir Shahmoradi, Wednesday 11:01 PM, Jan 15, 2014, ICMB, UT Austin
# Amir Shahmoradi, Friday 8:27 PM, Nov 8, 2013

#INPUT FILES:

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

entropy_HP = read.table('entropies/HP.entropy', header=T)
map_1RD8_AB = read.table('molecular_dynamics/map/HP_1RD8_AB.map',header=T)
entropy_1RD8_AB = entropy_HP$entropy[!is.na(map_1RD8_AB$pdb_pos)]; entropy_1RD8_AB = as.data.frame(entropy_1RD8_AB)
rsa_1RD8_AB = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1RD8_AB_sum.rsa',header=T)
rmsf_1RD8_AB = read.table('molecular_dynamics/Amber/postproc/rmsf/1RD8_AB_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_1RD8_AB = read.table('molecular_dynamics/Amber/postproc/dihedrals/1RD8_AB.dihedral',header=T)
cn13_1RD8_AB = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1RD8_AB_sum.cn13',header=T)
wcn_1RD8_AB = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1RD8_AB_sum.wcn',header=T)
desent_1RD8_AB = read.table('protein_design/entropies/1RD8_AB_0.6.desent',header=T)
bfca_1RD8_AB = read.table('bfactors/1RD8_AB_CA.bfactor',header=T,comment.char="")
rmsfHS_1RD8_AB = read.table('homologous_structures/rmsf/1RD8_AB_HS.rmsf',header=T,comment.char="")
omega_HP = read.table('evolutionary_rates/siterates_REL/H1N1_HA.txt', header=F, col.names = c("res_num", "omega"))
omega_1RD8_AB = omega_HP$omega[!is.na(map_1RD8_AB$pdb_pos)]; omega_1RD8_AB = as.data.frame(omega_1RD8_AB)

data_1RD8_AB = data.frame(protein="1RD8_AB", res_num=dihedral_1RD8_AB$res_Num, res_name=dihedral_1RD8_AB$res_name,
                       omega=omega_1RD8_AB$omega, entropy=entropy_1RD8_AB$entropy_1RD8_AB, desent=desent_1RD8_AB$entropy,
                       rsa_cr=rsa_1RD8_AB$CRYSTAL_RSA, rsa_avg_md=rsa_1RD8_AB$MEAN_RSA, rsa_var_md=rsa_1RD8_AB$VAR_RSA, 
                       rmsf_avg_md=rmsf_1RD8_AB$AvgRMSD, rmsf_std_md=rmsf_1RD8_AB$Stdev, 
                       phi_var_md=dihedral_1RD8_AB$var_phi, psi_var_md=dihedral_1RD8_AB$var_psi, chi1_var_md=dihedral_1RD8_AB$var_chi1,
                       cn13_cr=cn13_1RD8_AB$CRYSTAL_CN, cn13_avg_md=cn13_1RD8_AB$MEAN_CN, cn13_var_md=cn13_1RD8_AB$VAR_CN,
                       wcn_cr=wcn_1RD8_AB$CRYSTAL_WCN, wcn_avg_md=wcn_1RD8_AB$MEAN_WCN, wcn_var_md=wcn_1RD8_AB$VAR_WCN,
					   bfca=bfca_1RD8_AB$bfactor, rmsfHS=rmsfHS_1RD8_AB$rmsf)
write.csv( data_1RD8_AB, "correlation_analysis/combined_data/data_1RD8_AB.csv", row.names=F )


entropy_WNPB = read.table('entropies/WNPB.entropy', header=T)
map_2FP7_B = read.table('molecular_dynamics/map/WNPB_2FP7_B.map',header=T)
entropy_2FP7_B = entropy_WNPB$entropy[!is.na(map_2FP7_B$pdb_pos)]; entropy_2FP7_B = as.data.frame(entropy_2FP7_B)
rsa_2FP7_B = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2FP7_B_sum.rsa',header=T)
rmsf_2FP7_B = read.table('molecular_dynamics/Amber/postproc/rmsf/2FP7_B_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_2FP7_B = read.table('molecular_dynamics/Amber/postproc/dihedrals/2FP7_B.dihedral',header=T)
cn13_2FP7_B = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2FP7_B_sum.cn13',header=T)
wcn_2FP7_B = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2FP7_B_sum.wcn',header=T)
desent_2FP7_B = read.table('protein_design/entropies/2FP7_B_0.6.desent',header=T)
bfca_2FP7_B = read.table('bfactors/2FP7_B_CA.bfactor',header=T,comment.char="")
rmsfHS_2FP7_B = read.table('homologous_structures/rmsf/2FP7_B_HS.rmsf',header=T,comment.char="")
omega_WNPB = read.table('evolutionary_rates/siterates_REL/westnile_chainb.txt', header=F, col.names = c("res_num", "omega"))
omega_2FP7_B = omega_WNPB$omega[!is.na(map_2FP7_B$pdb_pos)]; omega_2FP7_B = as.data.frame(omega_2FP7_B)

data_2FP7_B = data.frame(protein="2FP7_B", res_num=dihedral_2FP7_B$res_Num, res_name=dihedral_2FP7_B$res_name,
                       omega=omega_2FP7_B$omega, entropy=entropy_2FP7_B$entropy_2FP7_B, desent=desent_2FP7_B$entropy,
                       rsa_cr=rsa_2FP7_B$CRYSTAL_RSA, rsa_avg_md=rsa_2FP7_B$MEAN_RSA, rsa_var_md=rsa_2FP7_B$VAR_RSA, 
                       rmsf_avg_md=rmsf_2FP7_B$AvgRMSD, rmsf_std_md=rmsf_2FP7_B$Stdev, 
                       phi_var_md=dihedral_2FP7_B$var_phi, psi_var_md=dihedral_2FP7_B$var_psi, chi1_var_md=dihedral_2FP7_B$var_chi1,
                       cn13_cr=cn13_2FP7_B$CRYSTAL_CN, cn13_avg_md=cn13_2FP7_B$MEAN_CN, cn13_var_md=cn13_2FP7_B$VAR_CN,
                       wcn_cr=wcn_2FP7_B$CRYSTAL_WCN, wcn_avg_md=wcn_2FP7_B$MEAN_WCN, wcn_var_md=wcn_2FP7_B$VAR_WCN,
					   bfca=bfca_2FP7_B$bfactor, rmsfHS=rmsfHS_2FP7_B$rmsf)
write.csv( data_2FP7_B, "correlation_analysis/combined_data/data_2FP7_B.csv", row.names=F )


entropy_JEHN = read.table('entropies/JEHN.entropy', header=T)
map_2Z83_A = read.table('molecular_dynamics/map/JEHN_2Z83_A.map',header=T)
entropy_2Z83_A = entropy_JEHN$entropy[!is.na(map_2Z83_A$pdb_pos)]; entropy_2Z83_A = as.data.frame(entropy_2Z83_A)
rsa_2Z83_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2Z83_A_sum.rsa',header=T)
rmsf_2Z83_A = read.table('molecular_dynamics/Amber/postproc/rmsf/2Z83_A_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_2Z83_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/2Z83_A.dihedral',header=T)
cn13_2Z83_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2Z83_A_sum.cn13',header=T)
wcn_2Z83_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2Z83_A_sum.wcn',header=T)
desent_2Z83_A = read.table('protein_design/entropies/2Z83_A_0.6.desent',header=T)
bfca_2Z83_A = read.table('bfactors/2Z83_A_CA.bfactor',header=T,comment.char="")
rmsfHS_2Z83_A = read.table('homologous_structures/rmsf/2Z83_A_HS.rmsf',header=T,comment.char="")
omega_JEHN = read.table('evolutionary_rates/siterates_REL/JEV.txt', header=F, col.names = c("res_num", "omega"))
omega_2Z83_A = omega_JEHN$omega[!is.na(map_2Z83_A$pdb_pos)]; omega_2Z83_A = as.data.frame(omega_2Z83_A)

data_2Z83_A = data.frame(protein="2Z83_A", res_num=dihedral_2Z83_A$res_Num, res_name=dihedral_2Z83_A$res_name,
                       omega=omega_2Z83_A$omega, entropy=entropy_2Z83_A$entropy_2Z83_A, desent=desent_2Z83_A$entropy,
                       rsa_cr=rsa_2Z83_A$CRYSTAL_RSA, rsa_avg_md=rsa_2Z83_A$MEAN_RSA, rsa_var_md=rsa_2Z83_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_2Z83_A$AvgRMSD, rmsf_std_md=rmsf_2Z83_A$Stdev, 
                       phi_var_md=dihedral_2Z83_A$var_phi, psi_var_md=dihedral_2Z83_A$var_psi, chi1_var_md=dihedral_2Z83_A$var_chi1,
                       cn13_cr=cn13_2Z83_A$CRYSTAL_CN, cn13_avg_md=cn13_2Z83_A$MEAN_CN, cn13_var_md=cn13_2Z83_A$VAR_CN,
                       wcn_cr=wcn_2Z83_A$CRYSTAL_WCN, wcn_avg_md=wcn_2Z83_A$MEAN_WCN, wcn_var_md=wcn_2Z83_A$VAR_WCN,
					   bfca=bfca_2Z83_A$bfactor, rmsfHS=rmsfHS_2Z83_A$rmsf)
write.csv( data_2Z83_A, "correlation_analysis/combined_data/data_2Z83_A.csv", row.names=F )


entropy_DPH = read.table('entropies/DPH.entropy', header=T)
map_2JLY_A = read.table('molecular_dynamics/map/DPH_2JLY_A.map',header=T)
entropy_2JLY_A = entropy_DPH$entropy[!is.na(map_2JLY_A$pdb_pos)]; entropy_2JLY_A = as.data.frame(entropy_2JLY_A)
rsa_2JLY_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2JLY_A_sum.rsa',header=T)
rmsf_2JLY_A = read.table('molecular_dynamics/Amber/postproc/rmsf/2JLY_A_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_2JLY_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/2JLY_A.dihedral',header=T)
cn13_2JLY_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2JLY_A_sum.cn13',header=T)
wcn_2JLY_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2JLY_A_sum.wcn',header=T)
desent_2JLY_A = read.table('protein_design/entropies/2JLY_A_0.6.desent',header=T)
bfca_2JLY_A = read.table('bfactors/2JLY_A_CA.bfactor',header=T,comment.char="")
rmsfHS_2JLY_A = read.table('homologous_structures/rmsf/2JLY_A_HS.rmsf',header=T,comment.char="")
omega_DPH = read.table('evolutionary_rates/siterates_REL/dengue_ns3.txt', header=F, col.names = c("res_num", "omega"))
omega_2JLY_A = omega_DPH$omega[!is.na(map_2JLY_A$pdb_pos)]; omega_2JLY_A = as.data.frame(omega_2JLY_A)
### UPDATE (APRIL 9, 2014): THIS PART OF THE CODE IS NOT NEEDED ANYMORE, SINCE DARIA HAS CALCULATED AND MAPPED ALL RMSF VALUES FOR 2JLY_A
### For CS RMSF, take the alignments from 2Z83, map 2JLY_A to JEHN alignemnt:
### mapped_RMSFHS_2Z83_A = c(); counter = 1
### for ( i in 1:nrow(map_2Z83_A) )
### {
### 	if ( is.na(map_2Z83_A$pdb_pos[i]) )
### 	{
### 		#cat( "i: ", as.character(i), "\n" )
### 		mapped_RMSFHS_2Z83_A[i] = NA
### 	}
### 	else
### 	{
### 		#print (counter)
### 		mapped_RMSFHS_2Z83_A[i] = rmsfHS_2Z83_A$rmsf[counter]
### 		counter = counter + 1
### 	}
### }
### mapped_RMSFHS_2Z83_A = cbind(map_2Z83_A,mapped_RMSFHS_2Z83_A)
### map_JEHN_2JLY_A = read.table('molecular_dynamics/map/JEHN_2JLY_A.map',header=T)
### csdata_2JLY_A = data.frame( aln_pos_2JLY_A = map_JEHN_2JLY_A$align_pos,  pdb_pos_2JLY_A = map_JEHN_2JLY_A$pdb_pos)
### csdata_2JLY_A = cbind(csdata_2JLY_A,mapped_RMSFHS_2Z83_A)
### rmsfHS_2JLY_A = csdata_2JLY_A$mapped_RMSFHS_2Z83_A[ !is.na(csdata_2JLY_A$pdb_pos_2JLY_A) ]

data_2JLY_A = data.frame(protein="2JLY_A", res_num=dihedral_2JLY_A$res_Num, res_name=dihedral_2JLY_A$res_name,
                       omega=omega_2JLY_A$omega, entropy=entropy_2JLY_A$entropy_2JLY_A, desent=desent_2JLY_A$entropy,
                       rsa_cr=rsa_2JLY_A$CRYSTAL_RSA, rsa_avg_md=rsa_2JLY_A$MEAN_RSA, rsa_var_md=rsa_2JLY_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_2JLY_A$AvgRMSD, rmsf_std_md=rmsf_2JLY_A$Stdev, 
                       phi_var_md=dihedral_2JLY_A$var_phi, psi_var_md=dihedral_2JLY_A$var_psi, chi1_var_md=dihedral_2JLY_A$var_chi1,
                       cn13_cr=cn13_2JLY_A$CRYSTAL_CN, cn13_avg_md=cn13_2JLY_A$MEAN_CN, cn13_var_md=cn13_2JLY_A$VAR_CN,
                       wcn_cr=wcn_2JLY_A$CRYSTAL_WCN, wcn_avg_md=wcn_2JLY_A$MEAN_WCN, wcn_var_md=wcn_2JLY_A$VAR_WCN,
					   bfca=bfca_2JLY_A$bfactor, rmsfHS=rmsfHS_2JLY_A$rmsf)
write.csv( data_2JLY_A, "correlation_analysis/combined_data/data_2JLY_A.csv", row.names=F )


entropy_HCP = read.table('entropies/HCP.entropy', header=T)
map_3GOL_A = read.table('molecular_dynamics/map/HCP_3GOL_A.map',header=T)
entropy_3GOL_A = entropy_HCP$entropy[!is.na(map_3GOL_A$pdb_pos)]; entropy_3GOL_A = as.data.frame(entropy_3GOL_A)
rsa_3GOL_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/3GOL_A_sum.rsa',header=T)
rmsf_3GOL_A = read.table('molecular_dynamics/Amber/postproc/rmsf/3GOL_A_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_3GOL_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/3GOL_A.dihedral',header=T)
cn13_3GOL_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/3GOL_A_sum.cn13',header=T)
wcn_3GOL_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/3GOL_A_sum.wcn',header=T)
desent_3GOL_A = read.table('protein_design/entropies/3GOL_A_0.6.desent',header=T)
bfca_3GOL_A = read.table('bfactors/3GOL_A_CA.bfactor',header=T,comment.char="")
rmsfHS_3GOL_A = read.table('homologous_structures/rmsf/3GOL_A_HS.rmsf',header=T,comment.char="")
omega_HCP = read.table('evolutionary_rates/siterates_REL/HCV_NS5B.txt', header=F, col.names = c("res_num", "omega"))
omega_3GOL_A = omega_HCP$omega[!is.na(map_3GOL_A$pdb_pos)]; omega_3GOL_A = as.data.frame(omega_3GOL_A)

data_3GOL_A = data.frame(protein="3GOL_A", res_num=dihedral_3GOL_A$res_Num, res_name=dihedral_3GOL_A$res_name,
                       omega=omega_3GOL_A$omega, entropy=entropy_3GOL_A$entropy_3GOL_A, desent=desent_3GOL_A$entropy,
                       rsa_cr=rsa_3GOL_A$CRYSTAL_RSA, rsa_avg_md=rsa_3GOL_A$MEAN_RSA, rsa_var_md=rsa_3GOL_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_3GOL_A$AvgRMSD, rmsf_std_md=rmsf_3GOL_A$Stdev, 
                       phi_var_md=dihedral_3GOL_A$var_phi, psi_var_md=dihedral_3GOL_A$var_psi, chi1_var_md=dihedral_3GOL_A$var_chi1,
                       cn13_cr=cn13_3GOL_A$CRYSTAL_CN, cn13_avg_md=cn13_3GOL_A$MEAN_CN, cn13_var_md=cn13_3GOL_A$VAR_CN,
                       wcn_cr=wcn_3GOL_A$CRYSTAL_WCN, wcn_avg_md=wcn_3GOL_A$MEAN_WCN, wcn_var_md=wcn_3GOL_A$VAR_WCN,
					   bfca=bfca_3GOL_A$bfactor, rmsfHS=rmsfHS_3GOL_A$rmsf)
write.csv( data_3GOL_A, "correlation_analysis/combined_data/data_3GOL_A.csv", row.names=F )


entropy_RVFVNP = read.table('entropies/RVFVNP.entropy', header=T)
map_3LYF_A = read.table('molecular_dynamics/map/RVFVNP_3LYF_A.map',header=T)
entropy_3LYF_A = entropy_RVFVNP$entropy[!is.na(map_3LYF_A$pdb_pos)]; entropy_3LYF_A = as.data.frame(entropy_3LYF_A)
rsa_3LYF_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/3LYF_A_sum.rsa',header=T)
rmsf_3LYF_A = read.table('molecular_dynamics/Amber/postproc/rmsf/3LYF_A_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_3LYF_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/3LYF_A.dihedral',header=T)
cn13_3LYF_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/3LYF_A_sum.cn13',header=T)
wcn_3LYF_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/3LYF_A_sum.wcn',header=T)
desent_3LYF_A = read.table('protein_design/entropies/3LYF_A_0.6.desent',header=T)
bfca_3LYF_A = read.table('bfactors/3LYF_A_CA.bfactor',header=T,comment.char="")
rmsfHS_3LYF_A = read.table('homologous_structures/rmsf/3LYF_A_HS.rmsf',header=T,comment.char="")
omega_RVFVNP = read.table('evolutionary_rates/siterates_REL/riftvalley.txt', header=F, col.names = c("res_num", "omega"))
omega_3LYF_A = omega_RVFVNP$omega[!is.na(map_3LYF_A$pdb_pos)]; omega_3LYF_A = as.data.frame(omega_3LYF_A)

data_3LYF_A = data.frame(protein="3LYF_A", res_num=dihedral_3LYF_A$res_Num, res_name=dihedral_3LYF_A$res_name,
                       omega=omega_3LYF_A$omega, entropy=entropy_3LYF_A$entropy_3LYF_A, desent=desent_3LYF_A$entropy,
                       rsa_cr=rsa_3LYF_A$CRYSTAL_RSA, rsa_avg_md=rsa_3LYF_A$MEAN_RSA, rsa_var_md=rsa_3LYF_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_3LYF_A$AvgRMSD, rmsf_std_md=rmsf_3LYF_A$Stdev, 
                       phi_var_md=dihedral_3LYF_A$var_phi, psi_var_md=dihedral_3LYF_A$var_psi, chi1_var_md=dihedral_3LYF_A$var_chi1,
                       cn13_cr=cn13_3LYF_A$CRYSTAL_CN, cn13_avg_md=cn13_3LYF_A$MEAN_CN, cn13_var_md=cn13_3LYF_A$VAR_CN,
                       wcn_cr=wcn_3LYF_A$CRYSTAL_WCN, wcn_avg_md=wcn_3LYF_A$MEAN_WCN, wcn_var_md=wcn_3LYF_A$VAR_WCN,
					   bfca=bfca_3LYF_A$bfactor, rmsfHS=rmsfHS_3LYF_A$rmsf)
write.csv( data_3LYF_A, "correlation_analysis/combined_data/data_3LYF_A.csv", row.names=F )

					   
entropy_CCHFN = read.table('entropies/CCHFN.entropy', header=T)
map_4AQF_B = read.table('molecular_dynamics/map/CCHFN_4AQF_B.map',header=T)
entropy_4AQF_B = entropy_CCHFN$entropy[!is.na(map_4AQF_B$pdb_pos)]; entropy_4AQF_B = as.data.frame(entropy_4AQF_B)
rsa_4AQF_B = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/4AQF_B_sum.rsa',header=T)
rmsf_4AQF_B = read.table('molecular_dynamics/Amber/postproc/rmsf/4AQF_B_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_4AQF_B = read.table('molecular_dynamics/Amber/postproc/dihedrals/4AQF_B.dihedral',header=T)
cn13_4AQF_B = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/4AQF_B_sum.cn13',header=T)
wcn_4AQF_B = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/4AQF_B_sum.wcn',header=T)
desent_4AQF_B = read.table('protein_design/entropies/4AQF_B_0.6.desent',header=T)
bfca_4AQF_B = read.table('bfactors/4AQF_B_CA.bfactor',header=T,comment.char="")
omega_CCHFN = read.table('evolutionary_rates/siterates_REL/crimean_congo.Nucleoprotein.txt', header=F, col.names = c("res_num", "omega"))
omega_4AQF_B = omega_CCHFN$omega[!is.na(map_4AQF_B$pdb_pos)]; omega_4AQF_B = as.data.frame(omega_4AQF_B)

data_4AQF_B = data.frame(protein="4AQF_B", res_num=dihedral_4AQF_B$res_Num, res_name=dihedral_4AQF_B$res_name,
                       omega=omega_4AQF_B$omega, entropy=entropy_4AQF_B$entropy_4AQF_B, desent=desent_4AQF_B$entropy,
                       rsa_cr=rsa_4AQF_B$CRYSTAL_RSA, rsa_avg_md=rsa_4AQF_B$MEAN_RSA, rsa_var_md=rsa_4AQF_B$VAR_RSA, 
                       rmsf_avg_md=rmsf_4AQF_B$AvgRMSD, rmsf_std_md=rmsf_4AQF_B$Stdev, 
                       phi_var_md=dihedral_4AQF_B$var_phi, psi_var_md=dihedral_4AQF_B$var_psi, chi1_var_md=dihedral_4AQF_B$var_chi1,
                       cn13_cr=cn13_4AQF_B$CRYSTAL_CN, cn13_avg_md=cn13_4AQF_B$MEAN_CN, cn13_var_md=cn13_4AQF_B$VAR_CN,
                       wcn_cr=wcn_4AQF_B$CRYSTAL_WCN, wcn_avg_md=wcn_4AQF_B$MEAN_WCN, wcn_var_md=wcn_4AQF_B$VAR_WCN,
					   bfca=bfca_4AQF_B$bfactor, rmsfHS=NA)
write.csv( data_4AQF_B, "correlation_analysis/combined_data/data_4AQF_B.csv", row.names=F )

					   
entropy_MRNABD = read.table('entropies/MRNABD.entropy', header=T)
map_4GHA_A = read.table('molecular_dynamics/map/MRNABD_4GHA_A.map',header=T)
entropy_4GHA_A = entropy_MRNABD$entropy[!is.na(map_4GHA_A$pdb_pos)]; entropy_4GHA_A = as.data.frame(entropy_4GHA_A)
rsa_4GHA_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/4GHA_A_sum.rsa',header=T)
rmsf_4GHA_A = read.table('molecular_dynamics/Amber/postproc/rmsf/4GHA_A_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_4GHA_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/4GHA_A.dihedral',header=T)
cn13_4GHA_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/4GHA_A_sum.cn13',header=T)
wcn_4GHA_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/4GHA_A_sum.wcn',header=T)
desent_4GHA_A = read.table('protein_design/entropies/4GHA_A_0.6.desent',header=T)
bfca_4GHA_A = read.table('bfactors/4GHA_A_CA.bfactor',header=T,comment.char="")
omega_MRNABD = read.table('evolutionary_rates/siterates_REL/marburg.vp35.txt', header=F, col.names = c("res_num", "omega"))
omega_4GHA_A = omega_MRNABD$omega[!is.na(map_4GHA_A$pdb_pos)]; omega_4GHA_A = as.data.frame(omega_4GHA_A)

data_4GHA_A = data.frame(protein="4GHA_A", res_num=dihedral_4GHA_A$res_Num, res_name=dihedral_4GHA_A$res_name,
                       omega=omega_4GHA_A$omega, entropy=entropy_4GHA_A$entropy_4GHA_A, desent=desent_4GHA_A$entropy,
                       rsa_cr=rsa_4GHA_A$CRYSTAL_RSA, rsa_avg_md=rsa_4GHA_A$MEAN_RSA, rsa_var_md=rsa_4GHA_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_4GHA_A$AvgRMSD, rmsf_std_md=rmsf_4GHA_A$Stdev, 
                       phi_var_md=dihedral_4GHA_A$var_phi, psi_var_md=dihedral_4GHA_A$var_psi, chi1_var_md=dihedral_4GHA_A$var_chi1,
                       cn13_cr=cn13_4GHA_A$CRYSTAL_CN, cn13_avg_md=cn13_4GHA_A$MEAN_CN, cn13_var_md=cn13_4GHA_A$VAR_CN,
                       wcn_cr=wcn_4GHA_A$CRYSTAL_WCN, wcn_avg_md=wcn_4GHA_A$MEAN_WCN, wcn_var_md=wcn_4GHA_A$VAR_WCN,
					   bfca=bfca_4GHA_A$bfactor, rmsfHS=NA)
write.csv( data_4GHA_A, "correlation_analysis/combined_data/data_4GHA_A.csv", row.names=F )


entropy_INP = read.table('entropies/INP.entropy', header=T)
map_4IRY_A = read.table('molecular_dynamics/map/INP_4IRY_A.map',header=T)
entropy_4IRY_A = entropy_INP$entropy[!is.na(map_4IRY_A$pdb_pos)]; entropy_4IRY_A = as.data.frame(entropy_4IRY_A)
rsa_4IRY_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/4IRY_A_sum.rsa',header=T)
rmsf_4IRY_A = read.table('molecular_dynamics/Amber/postproc/rmsf/4IRY_A_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_4IRY_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/4IRY_A.dihedral',header=T)
cn13_4IRY_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/4IRY_A_sum.cn13',header=T)
wcn_4IRY_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/4IRY_A_sum.wcn',header=T)
desent_4IRY_A = read.table('protein_design/entropies/4IRY_A_0.6.desent',header=T)
bfca_4IRY_A = read.table('bfactors/4IRY_A_CA.bfactor',header=T,comment.char="")
omega_INP = read.table('evolutionary_rates/siterates_REL/H1N1_NP.txt', header=F, col.names = c("res_num", "omega"))
omega_4IRY_A = omega_INP$omega[!is.na(map_4IRY_A$pdb_pos)]; omega_4IRY_A = as.data.frame(omega_4IRY_A)

data_4IRY_A = data.frame(protein="4IRY_A", res_num=dihedral_4IRY_A$res_Num, res_name=dihedral_4IRY_A$res_name,
                       omega=omega_4IRY_A$omega, entropy=entropy_4IRY_A$entropy_4IRY_A, desent=desent_4IRY_A$entropy,
                       rsa_cr=rsa_4IRY_A$CRYSTAL_RSA, rsa_avg_md=rsa_4IRY_A$MEAN_RSA, rsa_var_md=rsa_4IRY_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_4IRY_A$AvgRMSD, rmsf_std_md=rmsf_4IRY_A$Stdev, 
                       phi_var_md=dihedral_4IRY_A$var_phi, psi_var_md=dihedral_4IRY_A$var_psi, chi1_var_md=dihedral_4IRY_A$var_chi1,
                       cn13_cr=cn13_4IRY_A$CRYSTAL_CN, cn13_avg_md=cn13_4IRY_A$MEAN_CN, cn13_var_md=cn13_4IRY_A$VAR_CN,
                       wcn_cr=wcn_4IRY_A$CRYSTAL_WCN, wcn_avg_md=wcn_4IRY_A$MEAN_WCN, wcn_var_md=wcn_4IRY_A$VAR_WCN,
					   bfca=bfca_4IRY_A$bfactor, rmsfHS=NA)
write.csv( data_4IRY_A, "correlation_analysis/combined_data/data_4IRY_A.csv", row.names=F )


# HOMOLOGOUS STRUCTURES (IN THE SAME FAMILY OF 3GOL_A.PDB)

entropy_HCP = read.table('entropies/HCP.entropy', header=T)
map_3GSZ_A = read.table('molecular_dynamics/map/HCP_3GSZ_A.map',header=T)
entropy_3GSZ_A = entropy_HCP$entropy[!is.na(map_3GSZ_A$pdb_pos)]; entropy_3GSZ_A = as.data.frame(entropy_3GSZ_A)
rsa_3GSZ_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/3GSZ_A_sum.rsa',header=T)
rmsf_3GSZ_A = read.table('molecular_dynamics/Amber/postproc/rmsf/3GSZ_A_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_3GSZ_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/3GSZ_A.dihedral',header=T)
cn13_3GSZ_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/3GSZ_A_sum.cn13',header=T)
wcn_3GSZ_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/3GSZ_A_sum.wcn',header=T)
desent_3GSZ_A = read.table('protein_design/entropies/3GSZ_A_0.6.desent',header=T)
bfca_3GSZ_A = read.table('bfactors/3GSZ_A_CA.bfactor',header=T,comment.char="")
omega_HCP = read.table('evolutionary_rates/siterates_REL/HCV_NS5B.txt', header=F, col.names = c("res_num", "omega"))
omega_3GSZ_A = omega_HCP$omega[!is.na(map_3GSZ_A$pdb_pos)]; omega_3GSZ_A = as.data.frame(omega_3GSZ_A)

data_3GSZ_A = data.frame(protein="3GSZ_A", res_num=dihedral_3GSZ_A$res_Num, res_name=dihedral_3GSZ_A$res_name,
                       omega=omega_3GSZ_A$omega, entropy=entropy_3GSZ_A$entropy_3GSZ_A, desent=desent_3GSZ_A$entropy,
                       rsa_cr=rsa_3GSZ_A$CRYSTAL_RSA, rsa_avg_md=rsa_3GSZ_A$MEAN_RSA, rsa_var_md=rsa_3GSZ_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_3GSZ_A$AvgRMSD, rmsf_std_md=rmsf_3GSZ_A$Stdev, 
                       phi_var_md=dihedral_3GSZ_A$var_phi, psi_var_md=dihedral_3GSZ_A$var_psi, chi1_var_md=dihedral_3GSZ_A$var_chi1,
                       cn13_cr=cn13_3GSZ_A$CRYSTAL_CN, cn13_avg_md=cn13_3GSZ_A$MEAN_CN, cn13_var_md=cn13_3GSZ_A$VAR_CN,
                       wcn_cr=wcn_3GSZ_A$CRYSTAL_WCN, wcn_avg_md=wcn_3GSZ_A$MEAN_WCN, wcn_var_md=wcn_3GSZ_A$VAR_WCN,
					   bfca=bfca_3GSZ_A$bfactor, rmsfHS=NA)
write.csv( data_3GSZ_A, "correlation_analysis/combined_data/data_3GSZ_A.csv", row.names=F )


entropy_HCP = read.table('entropies/HCP.entropy', header=T)
map_3I5K_A = read.table('molecular_dynamics/map/HCP_3I5K_A.map',header=T)
entropy_3I5K_A = entropy_HCP$entropy[!is.na(map_3I5K_A$pdb_pos)]; entropy_3I5K_A = as.data.frame(entropy_3I5K_A)
rsa_3I5K_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/3I5K_A_sum.rsa',header=T)
rmsf_3I5K_A = read.table('molecular_dynamics/Amber/postproc/rmsf/3I5K_A_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_3I5K_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/3I5K_A.dihedral',header=T)
cn13_3I5K_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/3I5K_A_sum.cn13',header=T)
wcn_3I5K_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/3I5K_A_sum.wcn',header=T)
desent_3I5K_A = read.table('protein_design/entropies/3I5K_A_0.6.desent',header=T)
bfca_3I5K_A = read.table('bfactors/3I5K_A_CA.bfactor',header=T,comment.char="")
omega_HCP = read.table('evolutionary_rates/siterates_REL/HCV_NS5B.txt', header=F, col.names = c("res_num", "omega"))
omega_3I5K_A = omega_HCP$omega[!is.na(map_3I5K_A$pdb_pos)]; omega_3I5K_A = as.data.frame(omega_3I5K_A)

data_3I5K_A = data.frame(protein="3I5K_A", res_num=dihedral_3I5K_A$res_Num, res_name=dihedral_3I5K_A$res_name,
                       omega=omega_3I5K_A$omega, entropy=entropy_3I5K_A$entropy_3I5K_A, desent=desent_3I5K_A$entropy,
                       rsa_cr=rsa_3I5K_A$CRYSTAL_RSA, rsa_avg_md=rsa_3I5K_A$MEAN_RSA, rsa_var_md=rsa_3I5K_A$VAR_RSA, 
                       rmsf_avg_md=rmsf_3I5K_A$AvgRMSD, rmsf_std_md=rmsf_3I5K_A$Stdev, 
                       phi_var_md=dihedral_3I5K_A$var_phi, psi_var_md=dihedral_3I5K_A$var_psi, chi1_var_md=dihedral_3I5K_A$var_chi1,
                       cn13_cr=cn13_3I5K_A$CRYSTAL_CN, cn13_avg_md=cn13_3I5K_A$MEAN_CN, cn13_var_md=cn13_3I5K_A$VAR_CN,
                       wcn_cr=wcn_3I5K_A$CRYSTAL_WCN, wcn_avg_md=wcn_3I5K_A$MEAN_WCN, wcn_var_md=wcn_3I5K_A$VAR_WCN,
					   bfca=bfca_3I5K_A$bfactor, rmsfHS=NA)
write.csv( data_3I5K_A, "correlation_analysis/combined_data/data_3I5K_A.csv", row.names=F )


# DIFFERENT TEMPERATURE

rsa_2JLY_A_temp_50 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2JLY_A_temp_50_sum.rsa',header=T)
rmsf_2JLY_A_temp_50 = read.table('molecular_dynamics/Amber/postproc/rmsf/2JLY_A_temp_50_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_2JLY_A_temp_50 = read.table('molecular_dynamics/Amber/postproc/dihedrals/2JLY_A_temp_50.dihedral',header=T)
cn13_2JLY_A_temp_50 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2JLY_A_temp_50_sum.cn13',header=T)
wcn_2JLY_A_temp_50 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2JLY_A_temp_50_sum.wcn',header=T)

data_2JLY_A_temp_50 = data.frame(protein="2JLY_A_temp_50", res_num=dihedral_2JLY_A$res_Num, res_name=dihedral_2JLY_A$res_name,
                       omega=omega_2JLY_A$omega, entropy=entropy_2JLY_A$entropy_2JLY_A, desent=desent_2JLY_A$entropy,
                       rsa_cr=rsa_2JLY_A_temp_50$CRYSTAL_RSA, rsa_avg_md=rsa_2JLY_A_temp_50$MEAN_RSA, rsa_var_md=rsa_2JLY_A_temp_50$VAR_RSA, 
                       rmsf_avg_md=rmsf_2JLY_A_temp_50$AvgRMSD, rmsf_std_md=rmsf_2JLY_A_temp_50$Stdev, 
                       phi_var_md=dihedral_2JLY_A_temp_50$var_phi, psi_var_md=dihedral_2JLY_A_temp_50$var_psi, chi1_var_md=dihedral_2JLY_A_temp_50$var_chi1,
                       cn13_cr=cn13_2JLY_A_temp_50$CRYSTAL_CN, cn13_avg_md=cn13_2JLY_A_temp_50$MEAN_CN, cn13_var_md=cn13_2JLY_A_temp_50$VAR_CN,
                       wcn_cr=wcn_2JLY_A_temp_50$CRYSTAL_WCN, wcn_avg_md=wcn_2JLY_A_temp_50$MEAN_WCN, wcn_var_md=wcn_2JLY_A_temp_50$VAR_WCN,
					   bfca=bfca_2JLY_A$bfactor, rmsfHS=NA)
write.csv( data_2JLY_A_temp_50, "correlation_analysis/combined_data/data_2JLY_A_temp_50.csv", row.names=F )


rsa_2JLY_A_temp_100 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2JLY_A_temp_100_sum.rsa',header=T)
rmsf_2JLY_A_temp_100 = read.table('molecular_dynamics/Amber/postproc/rmsf/2JLY_A_temp_100_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_2JLY_A_temp_100 = read.table('molecular_dynamics/Amber/postproc/dihedrals/2JLY_A_temp_100.dihedral',header=T)
cn13_2JLY_A_temp_100 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2JLY_A_temp_100_sum.cn13',header=T)
wcn_2JLY_A_temp_100 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2JLY_A_temp_100_sum.wcn',header=T)

data_2JLY_A_temp_100 = data.frame(protein="2JLY_A_temp_100", res_num=dihedral_2JLY_A$res_Num, res_name=dihedral_2JLY_A$res_name,
                       omega=omega_2JLY_A$omega, entropy=entropy_2JLY_A$entropy_2JLY_A, desent=desent_2JLY_A$entropy,
                       rsa_cr=rsa_2JLY_A_temp_100$CRYSTAL_RSA, rsa_avg_md=rsa_2JLY_A_temp_100$MEAN_RSA, rsa_var_md=rsa_2JLY_A_temp_100$VAR_RSA, 
                       rmsf_avg_md=rmsf_2JLY_A_temp_100$AvgRMSD, rmsf_std_md=rmsf_2JLY_A_temp_100$Stdev, 
                       phi_var_md=dihedral_2JLY_A_temp_100$var_phi, psi_var_md=dihedral_2JLY_A_temp_100$var_psi, chi1_var_md=dihedral_2JLY_A_temp_100$var_chi1,
                       cn13_cr=cn13_2JLY_A_temp_100$CRYSTAL_CN, cn13_avg_md=cn13_2JLY_A_temp_100$MEAN_CN, cn13_var_md=cn13_2JLY_A_temp_100$VAR_CN,
                       wcn_cr=wcn_2JLY_A_temp_100$CRYSTAL_WCN, wcn_avg_md=wcn_2JLY_A_temp_100$MEAN_WCN, wcn_var_md=wcn_2JLY_A_temp_100$VAR_WCN,
					   bfca=bfca_2JLY_A$bfactor, rmsfHS=NA)
write.csv( data_2JLY_A_temp_100, "correlation_analysis/combined_data/data_2JLY_A_temp_100.csv", row.names=F )


rsa_2JLY_A_temp_200 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2JLY_A_temp_200_sum.rsa',header=T)
rmsf_2JLY_A_temp_200 = read.table('molecular_dynamics/Amber/postproc/rmsf/2JLY_A_temp_200_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_2JLY_A_temp_200 = read.table('molecular_dynamics/Amber/postproc/dihedrals/2JLY_A_temp_200.dihedral',header=T)
cn13_2JLY_A_temp_200 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2JLY_A_temp_200_sum.cn13',header=T)
wcn_2JLY_A_temp_200 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2JLY_A_temp_200_sum.wcn',header=T)

data_2JLY_A_temp_200 = data.frame(protein="2JLY_A_temp_200", res_num=dihedral_2JLY_A$res_Num, res_name=dihedral_2JLY_A$res_name,
                       omega=omega_2JLY_A$omega, entropy=entropy_2JLY_A$entropy_2JLY_A, desent=desent_2JLY_A$entropy,
                       rsa_cr=rsa_2JLY_A_temp_200$CRYSTAL_RSA, rsa_avg_md=rsa_2JLY_A_temp_200$MEAN_RSA, rsa_var_md=rsa_2JLY_A_temp_200$VAR_RSA, 
                       rmsf_avg_md=rmsf_2JLY_A_temp_200$AvgRMSD, rmsf_std_md=rmsf_2JLY_A_temp_200$Stdev, 
                       phi_var_md=dihedral_2JLY_A_temp_200$var_phi, psi_var_md=dihedral_2JLY_A_temp_200$var_psi, chi1_var_md=dihedral_2JLY_A_temp_200$var_chi1,
                       cn13_cr=cn13_2JLY_A_temp_200$CRYSTAL_CN, cn13_avg_md=cn13_2JLY_A_temp_200$MEAN_CN, cn13_var_md=cn13_2JLY_A_temp_200$VAR_CN,
                       wcn_cr=wcn_2JLY_A_temp_200$CRYSTAL_WCN, wcn_avg_md=wcn_2JLY_A_temp_200$MEAN_WCN, wcn_var_md=wcn_2JLY_A_temp_200$VAR_WCN,
					   bfca=bfca_2JLY_A$bfactor, rmsfHS=NA)
write.csv( data_2JLY_A_temp_200, "correlation_analysis/combined_data/data_2JLY_A_temp_200.csv", row.names=F )


rsa_2JLY_A_temp_450 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/2JLY_A_temp_450_sum.rsa',header=T)
rmsf_2JLY_A_temp_450 = read.table('molecular_dynamics/Amber/postproc/rmsf/2JLY_A_temp_450_Cpdb_CA.rmsf',header=T,comment.char="")
dihedral_2JLY_A_temp_450 = read.table('molecular_dynamics/Amber/postproc/dihedrals/2JLY_A_temp_450.dihedral',header=T)
cn13_2JLY_A_temp_450 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/2JLY_A_temp_450_sum.cn13',header=T)
wcn_2JLY_A_temp_450 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/2JLY_A_temp_450_sum.wcn',header=T)

data_2JLY_A_temp_450 = data.frame(protein="2JLY_A_temp_450", res_num=dihedral_2JLY_A$res_Num, res_name=dihedral_2JLY_A$res_name,
                       omega=omega_2JLY_A$omega, entropy=entropy_2JLY_A$entropy_2JLY_A, desent=desent_2JLY_A$entropy,
                       rsa_cr=rsa_2JLY_A_temp_450$CRYSTAL_RSA, rsa_avg_md=rsa_2JLY_A_temp_450$MEAN_RSA, rsa_var_md=rsa_2JLY_A_temp_450$VAR_RSA, 
                       rmsf_avg_md=rmsf_2JLY_A_temp_450$AvgRMSD, rmsf_std_md=rmsf_2JLY_A_temp_450$Stdev, 
                       phi_var_md=dihedral_2JLY_A_temp_450$var_phi, psi_var_md=dihedral_2JLY_A_temp_450$var_psi, chi1_var_md=dihedral_2JLY_A_temp_450$var_chi1,
                       cn13_cr=cn13_2JLY_A_temp_450$CRYSTAL_CN, cn13_avg_md=cn13_2JLY_A_temp_450$MEAN_CN, cn13_var_md=cn13_2JLY_A_temp_450$VAR_CN,
                       wcn_cr=wcn_2JLY_A_temp_450$CRYSTAL_WCN, wcn_avg_md=wcn_2JLY_A_temp_450$MEAN_WCN, wcn_var_md=wcn_2JLY_A_temp_450$VAR_WCN,
					   bfca=bfca_2JLY_A$bfactor, rmsfHS=NA)
write.csv( data_2JLY_A_temp_450, "correlation_analysis/combined_data/data_2JLY_A_temp_450.csv", row.names=F )


# OTHER

#rmsf_1AOR_A_temp_373 = read.table('molecular_dynamics/Amber/postproc/rmsf/1AOR_A_temp_373_rmsf_ref_inpcrd.txt',header=T,comment.char="")
#rsa_1AOR_A_temp_373 = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1AOR_A_temp_373_sum.rsa',header=T)
#dihedral_1AOR_A_temp_373 = read.table('molecular_dynamics/Amber/postproc/dihedrals/1AOR_A_temp_373_dihedrals.txt',header=T)
#cn13_1AOR_A_temp_373 = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1AOR_A_temp_373_sum.cn13',header=T)
#wcn_1AOR_A_temp_373 = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1AOR_A_temp_373_sum.wcn',header=T)
#
#data_1AOR_A_temp_373 = data.frame(protein="1AOR_A", res_num=dihedral_1AOR_A$res_Num, res_name=dihedral_1AOR_A$res_name,
#                       rsa_cr=rsa_1AOR_A_temp_373$CRYSTAL_RSA, rsa_avg_md=rsa_1AOR_A_temp_373$MEAN_RSA, rsa_var_md=rsa_1AOR_A_temp_373$VAR_RSA, 
#                       rmsf_avg_md=rmsf_1AOR_A_temp_373$AvgRMSD, rmsf_std_md=rmsf_1AOR_A_temp_373$Stdev, 
#                       phi_var_md=dihedral_1AOR_A_temp_373$var_phi, psi_var_md=dihedral_1AOR_A_temp_373$var_psi, chi1_var_md=dihedral_1AOR_A_temp_373$var_chi1,
#                       cn13_cr=cn13_1AOR_A_temp_373$CRYSTAL_CN, cn13_avg_md=cn13_1AOR_A_temp_373$MEAN_CN, cn13_var_md=cn13_1AOR_A_temp_373$VAR_CN,
#                       wcn_cr=wcn_1AOR_A_temp_373$CRYSTAL_WCN, wcn_avg_md=wcn_1AOR_A_temp_373$MEAN_WCN, wcn_var_md=wcn_1AOR_A_temp_373$VAR_WCN)

# THERMOPHILIC STRUCTURES

#rmsf_1AJ8_A = read.table('molecular_dynamics/Amber/postproc/rmsf/1AJ8_A_rmsf_ref_inpcrd.txt',header=T,comment.char="")
#rsa_1AJ8_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1AJ8_A_sum.rsa',header=T)
#dihedral_1AJ8_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/1AJ8_A_dihedrals.txt',header=T)
#cn13_1AJ8_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1AJ8_A_sum.cn13',header=T)
#wcn_1AJ8_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1AJ8_A_sum.wcn',header=T)
#
#data_1AJ8_A = data.frame(protein="1AJ8_A", res_num=dihedral_1AJ8_A$res_Num, res_name=dihedral_1AJ8_A$res_name,
#                       rsa_cr=rsa_1AJ8_A$CRYSTAL_RSA, rsa_avg_md=rsa_1AJ8_A$MEAN_RSA, rsa_var_md=rsa_1AJ8_A$VAR_RSA, 
#                       rmsf_avg_md=rmsf_1AJ8_A$AvgRMSD, rmsf_std_md=rmsf_1AJ8_A$Stdev, 
#                       phi_var_md=dihedral_1AJ8_A$var_phi, psi_var_md=dihedral_1AJ8_A$var_psi, chi1_var_md=dihedral_1AJ8_A$var_chi1,
#                       cn13_cr=cn13_1AJ8_A$CRYSTAL_CN, cn13_avg_md=cn13_1AJ8_A$MEAN_CN, cn13_var_md=cn13_1AJ8_A$VAR_CN,
#                       wcn_cr=wcn_1AJ8_A$CRYSTAL_WCN, wcn_avg_md=wcn_1AJ8_A$MEAN_WCN, wcn_var_md=wcn_1AJ8_A$VAR_WCN)
#
#rmsf_1AOR_A = read.table('molecular_dynamics/Amber/postproc/rmsf/1AOR_A_rmsf_ref_inpcrd.txt',header=T,comment.char="")
#rsa_1AOR_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1AOR_A_sum.rsa',header=T)
#dihedral_1AOR_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/1AOR_A_dihedrals.txt',header=T)
#cn13_1AOR_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1AOR_A_sum.cn13',header=T)
#wcn_1AOR_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1AOR_A_sum.wcn',header=T)
#
#data_1AOR_A = data.frame(protein="1AOR_A", res_num=dihedral_1AOR_A$res_Num, res_name=dihedral_1AOR_A$res_name,
#                       rsa_cr=rsa_1AOR_A$CRYSTAL_RSA, rsa_avg_md=rsa_1AOR_A$MEAN_RSA, rsa_var_md=rsa_1AOR_A$VAR_RSA, 
#                       rmsf_avg_md=rmsf_1AOR_A$AvgRMSD, rmsf_std_md=rmsf_1AOR_A$Stdev, 
#                       phi_var_md=dihedral_1AOR_A$var_phi, psi_var_md=dihedral_1AOR_A$var_psi, chi1_var_md=dihedral_1AOR_A$var_chi1,
#                       cn13_cr=cn13_1AOR_A$CRYSTAL_CN, cn13_avg_md=cn13_1AOR_A$MEAN_CN, cn13_var_md=cn13_1AOR_A$VAR_CN,
#                       wcn_cr=wcn_1AOR_A$CRYSTAL_WCN, wcn_avg_md=wcn_1AOR_A$MEAN_WCN, wcn_var_md=wcn_1AOR_A$VAR_WCN)
#
#rmsf_1CTS_A = read.table('molecular_dynamics/Amber/postproc/rmsf/1CTS_A_rmsf_ref_inpcrd.txt',header=T,comment.char="")
#rsa_1CTS_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1CTS_A_sum.rsa',header=T)
#dihedral_1CTS_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/1CTS_A_dihedrals.txt',header=T)
#cn13_1CTS_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1CTS_A_sum.cn13',header=T)
#wcn_1CTS_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1CTS_A_sum.wcn',header=T)
#
#data_1CTS_A = data.frame(protein="1CTS_A", res_num=dihedral_1CTS_A$res_Num, res_name=dihedral_1CTS_A$res_name,
#                       rsa_cr=rsa_1CTS_A$CRYSTAL_RSA, rsa_avg_md=rsa_1CTS_A$MEAN_RSA, rsa_var_md=rsa_1CTS_A$VAR_RSA, 
#                       rmsf_avg_md=rmsf_1CTS_A$AvgRMSD, rmsf_std_md=rmsf_1CTS_A$Stdev, 
#                       phi_var_md=dihedral_1CTS_A$var_phi, psi_var_md=dihedral_1CTS_A$var_psi, chi1_var_md=dihedral_1CTS_A$var_chi1,
#                       cn13_cr=cn13_1CTS_A$CRYSTAL_CN, cn13_avg_md=cn13_1CTS_A$MEAN_CN, cn13_var_md=cn13_1CTS_A$VAR_CN,
#                       wcn_cr=wcn_1CTS_A$CRYSTAL_WCN, wcn_avg_md=wcn_1CTS_A$MEAN_WCN, wcn_var_md=wcn_1CTS_A$VAR_WCN)
#
#rmsf_1MP9_A = read.table('molecular_dynamics/Amber/postproc/rmsf/1MP9_A_rmsf_ref_inpcrd.txt',header=T,comment.char="")
#rsa_1MP9_A = read.table('molecular_dynamics/Amber/postproc/rsa/summaries/1MP9_A_sum.rsa',header=T)
#dihedral_1MP9_A = read.table('molecular_dynamics/Amber/postproc/dihedrals/1MP9_A_dihedrals.txt',header=T)
#cn13_1MP9_A = read.table('molecular_dynamics/Amber/postproc/cn13/summaries/1MP9_A_sum.cn13',header=T)
#wcn_1MP9_A = read.table('molecular_dynamics/Amber/postproc/wcn/summaries/1MP9_A_sum.wcn',header=T)
#
#data_1MP9_A = data.frame(protein="1MP9_A", res_num=dihedral_1MP9_A$res_Num, res_name=dihedral_1MP9_A$res_name,
#                       rsa_cr=rsa_1MP9_A$CRYSTAL_RSA, rsa_avg_md=rsa_1MP9_A$MEAN_RSA, rsa_var_md=rsa_1MP9_A$VAR_RSA, 
#                       rmsf_avg_md=rmsf_1MP9_A$AvgRMSD, rmsf_std_md=rmsf_1MP9_A$Stdev, 
#                       phi_var_md=dihedral_1MP9_A$var_phi, psi_var_md=dihedral_1MP9_A$var_psi, chi1_var_md=dihedral_1MP9_A$var_chi1,
#                       cn13_cr=cn13_1MP9_A$CRYSTAL_CN, cn13_avg_md=cn13_1MP9_A$MEAN_CN, cn13_var_md=cn13_1MP9_A$VAR_CN,
#                       wcn_cr=wcn_1MP9_A$CRYSTAL_WCN, wcn_avg_md=wcn_1MP9_A$MEAN_WCN, wcn_var_md=wcn_1MP9_A$VAR_WCN)
#

