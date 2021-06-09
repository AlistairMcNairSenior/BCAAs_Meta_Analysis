
# Clean up the R environment
rm(list=ls())

# Set the WD
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/2018 BCAA Meta-analysis/Analysis" # Work iMac
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/2018 BCAA Meta-analysis/Analysis" # Work Macbook
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/2018 BCAA Meta-analysis/Analysis" # Home iMac
setwd(paste0(wd, "/Analyses"))

library(metafor)
library(plyr)
library(corpcor)
library(ape)
library(splines)
source("0.Headers.R")

# Load the results files, and format in to 2 part lists to make the tables
load("traits_pt1.Rdata")
traits_1<-traits
load("traits_pts2.3.Rdata")
traits_2<-traits
traits<-list(traits_1, traits_2)

load("data_list_pt1.Rdata")
data_1<-data_list
load("data_list_pts2.3.Rdata")
data_2<-data_list
data_list<-list(data_1, data_2)

load("TSV_list_pt1.Rdata")
TSV_1<-TSV_list
load("TSV_list_pts2.3.Rdata")
TSV_2<-TSV_list
TSV_list<-list(TSV_1, TSV_2)

load("MA_pool_list_pt1.Rdata")
MA_1<-MA_pool_list
load("MA_pool_list_pts2.3.Rdata")
MA_2<-MA_pool_list
MA_pool_list<-list(MA_1, MA_2)

load("ER_pool_list_pt1.Rdata")
ER_1<-ER_pool_list
load("ER_pool_list_pts2.3.Rdata")
ER_2<-ER_pool_list
ER_pool_list<-list(ER_1, ER_2)

load("TF_pool_list_pt1.Rdata")
TF_1<-TF_pool_list
load("TF_pool_list_pts2.3.Rdata")
TF_2<-TF_pool_list
TF_pool_list<-list(TF_1, TF_2)

load("MR_nutri_list_pt1.Rdata")
MR_nut_1<-MR_nutri_pool_list
load("MR_nutri_list_pts2.3.Rdata")
MR_nut_2<-MR_nutri_pool_list
MR_nutri_pool_list<-list(MR_nut_1, MR_nut_2)

load("MR_other_pool_list_pt1.Rdata")
MR_other_1<-MR_other_pool_list
load("MR_other_pool_list_pts2.3.Rdata")
MR_other_2<-MR_other_pool_list
MR_other_pool_list<-list(MR_other_1, MR_other_2)

# Change the WD to the tables folder
setwd(paste0(wd, "/tables"))

# Summary Table
for(p in 1:length(traits)){
	for(t in 1:length(traits[[p]])){
		
		# Grab trait t from part p 
		data<-data_list[[p]][[t]][[1]]
		
		# Generate a summary
		summary_t<-data.frame(trait=traits[[p]][[t]], n_articles=length(unique(data$Art_ID)), n_experiments=length(unique(data$Experimental_Unit)), n_diets=length(unique(c(data$Group, data$exp_Group))), n_lnRR=dim(data)[1], prop_mouse=mean(data$Species == "M"))
		
		# Combine the the results
		if(p == 1 & t==1){
			summary<-summary_t
		}else{
			summary<-rbind(summary, summary_t)
		}
	}
}

# Save the summary table
write.table(summary, file="summary_table.csv", sep=",", row.names=F, col.names=names(summary))

# Meta-Analysis Results
for(p in 1:length(traits)){
	for(t in 1:length(traits[[p]])){
		
		# Grab trait t from part p 
		model<-MA_pool_list[[p]][[t]][[1]]
		TSV<-TSV_list[[p]][[t]][[1]]
		tau2<-sum(model[[2]][,c(1,2)])
		I2_total<-tau2/(tau2+TSV) * 100
		I2_Exp<-model[[2]][,1]/(tau2+TSV) * 100
		
		# Generate a summary
		summary_t<-data.frame(trait=traits[[p]][[t]], Est.=model[[1]][,1], LCL=model[[1]][,1]-1.96*sqrt(model[[1]][,2]), UCL=model[[1]][,1]+1.96*sqrt(model[[1]][,2]), sigma2_Exp=model[[2]][,1], sigma2_Resid=model[[2]][,2], I2_total=I2_total, I2_Exp=I2_Exp)
		
		# Combine the the results
		if(p == 1 & t==1){
			summary<-summary_t
		}else{
			summary<-rbind(summary, summary_t)
		}
	}
}

# Save the summary table
write.table(summary, file="MA_table.csv", sep=",", row.names=F, col.names=names(summary))

# PB Results
for(p in 1:length(traits)){
	for(t in 1:length(traits[[p]])){
		
		# Grab trait t from part p 
		ER_model<-ER_pool_list[[p]][[t]][[1]]
		TF_model<-TF_pool_list[[p]][[t]][[1]]
		
		# Results table
		summary_t<-data.frame(trait=traits[[p]][[t]], ER_coef=ER_model[[1]][2,1], ER_LCL=ER_model[[1]][2,1] - 1.96*sqrt(ER_model[[1]][2,3]), ER_UCL=ER_model[[1]][2,1] + 1.96*sqrt(ER_model[[1]][2,3]), n_missing=TF_model[[3]][1], prop_right=TF_model[[3]][2], TF_coef=TF_model[[1]][1], TF_LCL=TF_model[[1]][1] - 1.96*sqrt(TF_model[[1]][2]), TF_UCL=TF_model[[1]][1] + 1.96*sqrt(TF_model[[1]][2]))
		
		# Combine the the results
		if(p == 1 & t==1){
			summary<-summary_t
		}else{
			summary<-rbind(summary, summary_t)
		}
	}
}

# Save the summary table
write.table(summary, file="PB_table.csv", sep=",", row.names=F, col.names=names(summary))

# 'Other' Meta-Regression Results
for(p in 1:length(traits)){
	for(t in 1:length(traits[[p]])){
		
		# Grab trait t from part p 
		models<-MR_other_pool_list[[p]][[t]][[1]]
		
		# Results table
		summary_t<-data.frame(trait=NA, Coef=NA, Est.=NA, LCL=NA, UCL=NA, Q_Mod=NA, Q_p=NA, m.models=NA)
		
		# Go through each model and pull out the results
		for(r in 1:length(models)){
			
			# pull out the rth_model
			model_r<-models[[r]]
			
			# Get the numebr of coefficients
			n_coef<-length(row.names(model_r[[1]]))
			
			# Check out that the model actually fitted
			if(is.na(model_r[[1]][1,1]) == F){
				
				# Get the SE
				se<-sqrt(diag(model_r[[1]][,-1]))
				
				# combine the results
				summary_t_r<-data.frame(trait=c(traits[[p]][[t]], rep("", n_coef-1)), Coef=row.names(model_r[[1]]), Est.=model_r[[1]][,1], LCL=model_r[[1]][,1]-1.96*se, UCL=model_r[[1]][,1]+1.96*se, Q_Mod=c(model_r[[2]][,3], rep("", n_coef-1)), Q_p=c(model_r[[2]][,4], rep("", n_coef-1)), m.models=model_r[[4]])	
				
				# Add to the other results for trait t
				summary_t<-rbind(summary_t, summary_t_r)	
			}
		}
		
		summary_t<-summary_t[-1,]

		# Combine the the results
		if(p == 1 & t==1){
			summary<-summary_t
		}else{
			summary<-rbind(summary, summary_t)
		}
	}
}

# Save the summary table
write.table(summary, file="MR_table.csv", sep=",", row.names=F, col.names=names(summary))

# AIC of nutritional models
for(p in 1:length(traits)){
	for(t in 1:length(traits[[p]])){
	
		# Get the AIC of the meta-analysis
		summary_t<-data.frame(trait=traits[[p]][[t]])	
		summary_t$Model<-0
		summary_t$Model_name<-"MA"
		summary_t$df<-dim(MA_pool_list[[p]][[t]][[1]][[1]])[1]	
		summary_t$AIC<-MA_pool_list[[p]][[t]][[1]][[3]]
		
		# Get AIC from MRs
		for(r in 1:length(MR_nutri_pool_list[[p]][[t]][[1]])){
			model_r<-MR_nutri_pool_list[[p]][[t]][[1]][[r]]
			if(row.names(model_r[[1]])[1] != "NoModel"){
				summary_t<-rbind(summary_t, data.frame(trait="", Model=r, Model_name=names(MR_nutri_pool_list[[p]][[t]][[1]])[r], df=length(row.names(model_r[[1]])), AIC=model_r[[3]]))
			}
		}
		
		# Sort on AIC and pick within 2 points the model with lowest DF and then AIC
		sort_summary_t<-summary_t[order(summary_t$AIC),]
		summary_t[,-1]<-sort_summary_t[,-1]
		summary_t$delta_AIC<-summary_t$AIC - summary_t$AIC[1]
		summary_t$selection<-0
		options<-summary_t[which(summary_t$delta_AIC < 2),]
		choice<-options$Model[which(options$df == min(options$df))[1]]
		summary_t$selection[which(summary_t$Model == choice)]<-1
		
		# Add all the traits together
		if(p == 1 & t == 1){
			summary<-summary_t
		}else{
			summary<-rbind(summary, summary_t)	
		}	
	}
}

write.table(summary, file="AIC_table.csv", sep=",", row.names=F, col.names=names(summary))


# AIC-favoured Nutritional Meta-Regression Results

AIC<-summary
choices<-AIC$Model[which(AIC$selection == 1)]
aic_models<-list()

counter<-1
for(p in 1:length(traits)){
	for(t in 1:length(traits[[p]])){
		
		# Grab trait t from part p 
		trait<-traits[[p]][[t]]
		
		# Get the model, and also save it
		if(choices[counter] != 0){
			model<-MR_nutri_pool_list[[p]][[t]][[1]][[choices[counter]]]
		}else{
			model<-MA_pool_list[[p]][[t]][[1]]
		}
		aic_models[[counter]]<-model
		names(aic_models)[counter]<-trait
		
		
		# Results table
		Ests<-model[[1]][,1]
		SE<-sqrt(diag(as.matrix(model[[1]][,-1])))
		
		summary_t<-data.frame(trait=c(trait, rep("", length(Ests)-1)), Coef=row.names(model[[1]]), Est.=Ests, LCL=Ests-1.96*SE, UCL=Ests+1.96*SE, Q_Mod=c(model[[2]][3], rep("", length(Ests)-1)), Q_p=c(model[[2]][4], rep("", length(Ests)-1)), m.models=c(model[[4]], rep("", length(Ests)-1)))

		# Combine the the results
		if(p == 1 & t==1){
			summary<-summary_t
		}else{
			summary<-rbind(summary, summary_t)
		}
		
		#Advan ce the counter
		counter<-counter+1
	}
}

# Save the summary table
write.table(summary, file="MR_nutri_table.csv", sep=",", row.names=F, col.names=names(summary))

# Save the AIC favoured models
setwd(paste0(wd, "/Analyses"))
save(aic_models, file="AIC_models.Rdata")
