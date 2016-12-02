source("profilingFunctions.R",local=TRUE)

# List of filenames to open
fnames_MDRMDP=list(
ura3_27="QFA0141_FitnessReport_DAL_ura3_27_SDM_rhk_CTGN_MDRMDP.txt",
lyp1_30="QFA0132_FitnessReport_DAL_lyp1_HLN_30_SDM_rhlk_CTGNH_MDRMDP.txt",
lyp1_33="QFA0132_FitnessReport_DAL_lyp1_HLN_33_SDM_rhlk_CTGNH_MDRMDP.txt",
yku70_375="QFA0139_FitnessReport_DAL_Yku70_37_5_SDM_rhk_CTGN_MDRMDP.txt",
cdc13_27="QFA0140_FitnessReport_DAL_cdc13-1_27_SDM_rhlk_CTGH_MDRMDP.txt",
stn1_33="QFA0136_FitnessReport_DAL_stn1-13_33_SDM_rhlk_CTGNH_MDRMDP.txt",
rfa3_30="QFA0131_FitnessReport_MJG_rfa3-313_30_SDM_rhlk_CTGNH_MDRMDP.txt",
cdc13_exo1_30="QFA0051_FitnessReport_DAL_cdc13-1_exo1D_30_SDM_rhlk_CTGNH_MDRMDP.txt",
cdc13_UD="QFA0140_FitnessReport_DAL_cdc13-1_UD_X3_SDM_rhlk_CTGH_MDRMDP.txt",
cdc13_rad9_UD="QFA0142_FitnessReport_APB_cdc13-1_rad9D_UD_X1_SDM_rhlk_CTGNH_MDRMDP.txt",
cdc13_rad9_27="QFA0142_FitnessReport_APB_cdc13-1_rad9D_27_SDM_rhlk_CTGNH_MDRMDP.txt"
)

fnames_nAUC=list(
ura3_27="QFA0141_FitnessReport_DAL_ura3_27_SDM_rhk_CTGN_nAUC.txt",
lyp1_30="QFA0132_FitnessReport_DAL_lyp1_HLN_30_SDM_rhlk_CTGNH_nAUC.txt",
lyp1_33="QFA0132_FitnessReport_DAL_lyp1_HLN_33_SDM_rhlk_CTGNH_nAUC.txt",
yku70_375="QFA0139_FitnessReport_DAL_Yku70_37_5_SDM_rhk_CTGN_nAUC.txt",
cdc13_27="QFA0140_FitnessReport_DAL_cdc13-1_27_SDM_rhlk_CTGH_nAUC.txt",
stn1_33="QFA0136_FitnessReport_DAL_stn1-13_33_SDM_rhlk_CTGNH_nAUC.txt",
rfa3_30="QFA0131_FitnessReport_MJG_rfa3-313_30_SDM_rhlk_CTGNH_nAUC.txt",
cdc13_exo1_30="QFA0051_FitnessReport_DAL_cdc13-1_exo1D_30_SDM_rhlk_CTGNH_nAUC.txt",
cdc13_UD="QFA0140_FitnessReport_DAL_cdc13-1_UD_X3_SDM_rhlk_CTGH_nAUC.txt",
cdc13_rad9_UD="QFA0142_FitnessReport_APB_cdc13-1 rad9D_UD23-36_8h_1x _SDM_rhlk_CTGH_nAUC.txt",
cdc13_rad9_27="QFA0142_FitnessReport_APB_cdc13-1 rad9D_27_SDM_rhlk_CTGH_nAUC.txt"
)

stripType="GeneStripped"
for (f in names(fnames_nAUC)){
  fnames_nAUC[[f]]=file.path(stripType,fnames_nAUC[[f]])
  fnames_MDRMDP[[f]]=file.path(stripType,fnames_MDRMDP[[f]])
}

flist=fnames_MDRMDP