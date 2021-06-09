
# Clean up the R environment
rm(list=ls())

# Set the WD
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/2018 BCAA Meta-analysis/Analysis" # Work iMac
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/2018 BCAA Meta-analysis/Analysis" # Work Macbook
setwd(paste0(wd, "/Analyses"))

# Load the libraries
library(mice)
library(metafor)
library(plyr)
library(corpcor)
library(ape)
library(splines)
library(tcltk)
library(doSNOW)
source("0.Headers.R")

# How many imputations will we be doing for missing data
n_imp<-20

################################################################
############################ PART 1 ############################
################################################################

# In Pt 1 we look at the effects of dietary BCAAs on circulating BCAA levels
traits_pt1<-c("circ_bcaa", "circ_iso", "circ_leu", "circ_val")

# List the traits to model
traits<-traits_pt1

# Meta-Regressions list
# For each of circulating amino-acid measure, we will compare that amino-acid (q), as well as P - selection by AIC
# Models will be fitted where n_effects >= params * 10

# 1) q in control: lin, spline. params = 2, df+1.
# 2) P in control: lin, spline. params = 2, df+1.
# 3) delta q: lin, spline. params = 2, df+1.
# 4) delta P: lin, spline. params = 2, df+1.
# 5) q in control * delta q: lin, spline. params = 4, df^2+1.
# 6) P in control * delta P: lin, spline. params = 4, df^2+1.

# In addition we look at seperate models Species, fasting, gavage
formulas_list_2<-list()
params_list_2<-list()
formulas_list_2[[1]]<-"Species"
formulas_list_2[[2]]<-"Time_fasted"
params_list_2[[1]]<-NA
params_list_2[[2]]<-NA


# The nutrients as listed in the dataset column names
qs<-c("BCAA", "Iso", "Leu", "Val")

# Creating formulas for all of these nutritional models
formulas_list_all<-list()
params_list<-list()

# Go through each and create a unique list for that amino acid measure
for(i in 1:length(qs)){
	
	# For each create a unique set of formulas
	formulas_list<-list()
	
	# Models 1 and 2: lin and spline version
	formulas_list[[1]]<-paste0("Diet_", qs[i], "_kJ_g")
	params_list[[1]]<-2
	formulas_list[[2]]<-paste0("bs(Diet_", qs[i], "_kJ_g, df=3)")
	params_list[[2]]<-4
	formulas_list[[3]]<-paste0("Diet_P_kJ_g")
	params_list[[3]]<-2
	formulas_list[[4]]<-paste0("bs(Diet_P_kJ_g, df=3)")
	params_list[[4]]<-4

	
	# Models 3 and 4: lin and spline 
	formulas_list[[5]]<-paste0("diff_", qs[i])
	params_list[[5]]<-2
	formulas_list[[6]]<-paste0("bs(diff_", qs[i], ", df=3)")
	params_list[[6]]<-4
	formulas_list[[7]]<-paste0("diff_P")
	params_list[[7]]<-2
	formulas_list[[8]]<-paste0("bs(diff_P, df=3)")
	params_list[[8]]<-4
	
	# Models 5 and 6: lin and spline
	formulas_list[[9]]<-paste0("Diet_", qs[i], "_kJ_g * diff_", qs[i])
	params_list[[9]]<-4
	formulas_list[[10]]<-paste0("bs(Diet_", qs[i], "_kJ_g, df=3) : bs(diff_", qs[i], ", df=3)")
	params_list[[10]]<-10
	formulas_list[[11]]<-paste0("Diet_P_kJ_g * diff_P")
	params_list[[11]]<-4
	formulas_list[[12]]<-paste0("bs(Diet_P_kJ_g, df=3) : bs(diff_P, df=3)")
	params_list[[12]]<-10
	
	# Save the ith set
	formulas_list_all[[i]]<-formulas_list
}


# I am going to run over 4 clusters (one per trait)
cl<-makeCluster(length(traits_pt1), outfile="")
registerDoSNOW(cl)

# run function for each part across the four cores
results<-foreach(p = 1:length(traits_pt1)) %dopar% {	
	output<-analyse_traits(traits=traits_pt1[p], nutri_formulas=formulas_list_all[[p]], nutri_params=params_list, other_formulas=formulas_list_2, other_params=params_list_2, n_imp=n_imp, wd=wd, core=p)
	return(output)
}


# Unparcel the results
data_list<-lapply(results, "[[", 1)
TSV_list<-lapply(results, "[[", 2)
MA_pool_list<-lapply(results, "[[", 3)
ER_pool_list<-lapply(results, "[[", 4)
TF_pool_list<-lapply(results, "[[", 5)
MR_nutri_pool_list<-lapply(results, "[[", 6)
MR_other_pool_list<-lapply(results, "[[", 7)

# Save
save(traits, file="traits_pt1.Rdata")
save(data_list, file="data_list_pt1.Rdata")
save(TSV_list, file="TSV_list_pt1.Rdata")
save(MA_pool_list, file="MA_pool_list_pt1.Rdata")
save(ER_pool_list, file="ER_pool_list_pt1.Rdata")
save(TF_pool_list, file="TF_pool_list_pt1.Rdata")
save(MR_nutri_pool_list, file="MR_nutri_list_pt1.Rdata")
save(MR_other_pool_list, file="MR_other_pool_list_pt1.Rdata")

################################################################
######################## PARTS 2 & 3 ###########################
################################################################

# In parts 2 and 3 we look at the effects of dietary BCAAs and other parameters on glucose homeostasis, body mass composition and food intake 

# Pt 2 - glucose homeostasis
traits_pt2<-c("glu_AUC", "ins_AUC", "plas_glu", "plas_ins", "HOMA")

# Pt 3 - body mass, compoisition and food intake
traits_pt3<-c("weight", "%_fat", "intake", "energy_intake")

# List of traits
traits<-c(traits_pt2, traits_pt3)

# Meta-Regressions list
# OK so we need to explore the follwing nutritional models - selection by AIC
# Models will be fitted where n_effects >= params * 10

# 1) BCAA in control: lin, spline. params = 2, df+1.
# 2) P in control: lin, spline. params = 2, df+1.
# 3) PC in control: lin, spline. params = 2, df+1.
# 4) delta BCAA: lin, spline. params = 2, df+1.
# 5) delta BCAA:non-BCAA: line, spine. params = 2, df+1.
# 6) BCAA in control * delta BCAA: lin, spline. params = 4, df^2+1.
# 7) P in control * delta BCAA: lin, spline. params = 4, df^2+1.
# 8) PC in control * delta BCAA: lin, spline. params = 4, df^2+1.
# 9) delta non_BCAA * delta BCAA: lin, spline. params = 4, df^2+1.

# In addition we look at seperate models Species, fasting, gavage

# Creating formulas for all of these models
formulas_list<-list()
params_list<-list()

# 1) BCAA in control: lin, spline. params = 2, df+1.
formulas_list[[1]]<-"Diet_BCAA_kJ_g"
formulas_list[[2]]<-"bs(Diet_BCAA_kJ_g, df=3)"
params_list[[1]]<-2
params_list[[2]]<-5

# 2) P in control: lin, spline. params = 2, df+1.
formulas_list[[3]]<-"Diet_P_kJ_g"
formulas_list[[4]]<-"bs(Diet_P_kJ_g, df=3)"
params_list[[3]]<-2
params_list[[4]]<-5

# 3) PC in control: lin, spline. params = 2, df+1.
formulas_list[[5]]<-"Diet_PC_kJ_g"
formulas_list[[6]]<-"bs(Diet_PC_kJ_g, df=3)"
params_list[[5]]<-2
params_list[[6]]<-5

# 4) delta_BCAA in control: lin, spline. params = 2, df+1.
formulas_list[[7]]<-"diff_BCAA"
formulas_list[[8]]<-"bs(diff_BCAA, df=3)"
params_list[[7]]<-2
params_list[[8]]<-5

# 4) delta_BCAA to non_BCAA: lin, spline. params = 2, df+1.
formulas_list[[9]]<-"diff_BCAA.Non_BCAA_kJ_g"
formulas_list[[10]]<-"bs(diff_BCAA.Non_BCAA_kJ_g, df=3)"
params_list[[9]]<-2
params_list[[10]]<-5

#6) BCAA in control * delta BCAA: lin, spline. params = 4, df^2+1.
formulas_list[[11]]<-"Diet_BCAA_kJ_g * diff_BCAA"
formulas_list[[12]]<-"bs(Diet_BCAA_kJ_g, df=3) : bs(diff_BCAA, df=3)"
params_list[[11]]<-4
params_list[[12]]<-17

#7) P in control * delta BCAA: lin, spline. params = 4, df^2+1.
formulas_list[[13]]<-"Diet_P_kJ_g * diff_BCAA"
formulas_list[[14]]<-"bs(Diet_P_kJ_g, df=3) : bs(diff_BCAA, df=3)"
params_list[[13]]<-4
params_list[[14]]<-17

#8) PC in control * delta BCAA: lin, spline. params = 4, df^2+1.
formulas_list[[15]]<-"Diet_PC_kJ_g * diff_BCAA"
formulas_list[[16]]<-"bs(Diet_PC_kJ_g, df=3) : bs(diff_BCAA, df=3)"
params_list[[15]]<-4
params_list[[16]]<-17

#9) delta non_BCAA * delta BCAA: lin, spline. params = 4, df^2+1.
formulas_list[[17]]<-"diff_nonBCAA * diff_BCAA"
formulas_list[[18]]<-"bs(diff_nonBCAA, df=3) : bs(diff_BCAA, df=3)"
params_list[[17]]<-4
params_list[[18]]<-17

# Other models
formulas_list_2<-list()
params_list_2<-list()
formulas_list_2[[1]]<-"Species"
# Note to self - I have done the models for time fasted, but for the traits in Pt3 these do not make any sense and so should be removed
formulas_list_2[[2]]<-"Time_fasted"
params_list_2[[1]]<-NA
params_list_2[[2]]<-NA

# I am going to run over 9 clusters (one per part)
cl<-makeCluster(length(traits), outfile="")
registerDoSNOW(cl)

# run function for each part across the four cores
results<-foreach(p = 1:length(traits)) %dopar% {
	output<-analyse_traits(traits=traits[p], nutri_formulas=formulas_list, nutri_params=params_list, other_formulas=formulas_list_2, other_params=params_list_2, n_imp=n_imp, wd=wd, core=p)
	return(output)
}

# Unparcel the results
data_list<-lapply(results, "[[", 1)
TSV_list<-lapply(results, "[[", 2)
MA_pool_list<-lapply(results, "[[", 3)
ER_pool_list<-lapply(results, "[[", 4)
TF_pool_list<-lapply(results, "[[", 5)
MR_nutri_pool_list<-lapply(results, "[[", 6)
MR_other_pool_list<-lapply(results, "[[", 7)

# Save
save(traits, file="traits_pts2.3.Rdata")
save(data_list, file="data_list_pts2.3.Rdata")
save(TSV_list, file="TSV_list_pts2.3.Rdata")
save(MA_pool_list, file="MA_pool_list_pts2.3.Rdata")
save(ER_pool_list, file="ER_pool_list_pts2.3.Rdata")
save(TF_pool_list, file="TF_pool_list_pts2.3.Rdata")
save(MR_nutri_pool_list, file="MR_nutri_list_pts2.3.Rdata")
save(MR_other_pool_list, file="MR_other_pool_list_pts2.3.Rdata")


