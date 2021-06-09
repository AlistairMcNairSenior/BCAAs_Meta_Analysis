
# Clean up the R environment
rm(list=ls())

# Set the WD
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/2018 BCAA Meta-analysis/Analysis" # Work iMac
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/2018 BCAA Meta-analysis/Analysis" # Work Macbook
#wd<-"/Users/asenior/Dropbox (Sydney Uni)/2018 BCAA Meta-analysis/Analysis" # Home iMac
setwd(paste0(wd, "/Analyses"))

library(metafor)
library(Cairo)
library(splines)
library(ggplot2)
library(sp)
library(metR)
library(gridExtra)
library(plyr)
library(corpcor)
library(ape)
source("0.Headers.R")

# Load the results files, and format in to 2 part lists to make the tables
load("traits_pt1.Rdata")
traits_1<-traits
load("traits_pts2.3.Rdata")
traits_2<-traits
# Note I am subsetting traits_2 for convenience to break apart parts 2 and 3 for plotting
traits<-list(traits_1, traits_2[c(1:5)], traits_2[c(6:9)])

load("data_list_pt1.Rdata")
data_1<-data_list
load("data_list_pts2.3.Rdata")
data_2<-data_list
# Note I am subsetting traits_2 for convenience to break apart parts 2 and 3 for plotting
data_list<-list(data_1, data_2[c(1:5)], data_2[c(6:9)])

load("TSV_list_pt1.Rdata")
TSV_1<-TSV_list
load("TSV_list_pts2.3.Rdata")
TSV_2<-TSV_list
TSV_list<-list(TSV_1, TSV_2[c(1:5)], TSV_2[c(6:9)])

load("MA_pool_list_pt1.Rdata")
MA_1<-MA_pool_list
load("MA_pool_list_pts2.3.Rdata")
MA_2<-MA_pool_list
MA_pool_list<-list(MA_1, MA_2[c(1:5)], MA_2[c(6:9)])

# load the AIC Favoured models
load("AIC_models.rdata")
aic_models
aic_models<-list(aic_models[c(1:4)], aic_models[c(5:9)], aic_models[c(10:13)])

# Load the Other MR models to make the fasting vs non-fasting plots
load("MR_other_pool_list_pt1.Rdata")
MR_1<-MR_other_pool_list
load("MR_other_pool_list_pts2.3.Rdata")
MR_2<-MR_other_pool_list

# We only want the fasting vs non-fasting results, and for circulating levels (pt 1) and also glucose metabolism traits (pt 2)
MR_fast_list<-list()
sub<-list()
for(i in 1:4){
	sub[[i]]<-MR_1[[i]][[1]]$Time_fasted
}
MR_fast_list[[1]]<-sub
for(i in 1:5){
	sub[[i]]<-MR_2[[i]][[1]]$Time_fasted
}
MR_fast_list[[2]]<-sub

# What are the types of models that we are looking at - need to know for plotting. Bubble plots for single dimensions, surfaces for 2 dimensions - if MA we just pass by		# BCAA			# ISO	# Leu		# Val	# GLU AUC	# Ins AUC	$PG	  # Plas Ins.  #HOMA			# Mass		# %F		# Food		# Energy	 
aic_types<-list(c("spline.1d", "lin.2d", "lin.2d", "lin.2d"), c("lin.2d", "spline.1d", "ma", "spline.1d", "spline.1d"), c("spline.1d", "lin.2d", "spline.1d", "spline.1d"))

# What was the df used in the splines
df<-3

# Set to the figures folder for output
setwd(paste0(wd, "/figures"))

# Lets make some nice titles for each trait
# traits
trait_titles<-list()
trait_titles[[1]]<-c("Plasma
BCAAs", "Plasma
Isoleucine", "Plasma
Leucine", "Plasma
Valine")
trait_titles[[2]]<-c("Glucose
AUC", "Insulin
AUC", "Plasma
Glucose", "Plasma
Insulin", "HOMA")
trait_titles[[3]]<-c("Body
Mass", "Percent
Fat Mass", "Food
Intake", "Energy
Intake")

# Make some nice names for the meta-regression moderators
pred_titles_pt1<-list()
pred_titles_pt1[[1]]<-"Control Diet BCAA (kJ/g)"
pred_titles_pt1[[2]]<-c("Control Diet Isoleucine (kJ/g)", "\U0394 Isoleucine (kJ/g)") 
pred_titles_pt1[[3]]<-c("Control Diet Leucine (kJ/g)", "\U0394 Leucine (kJ/g)")
pred_titles_pt1[[4]]<-c("Control Diet Valine (kJ/g)", "\U0394 Valine (kJ/g)") 

pred_titles_pt2<-list()
pred_titles_pt2[[1]]<-c("\U0394 Non-BCAA (kJ/g)", "\U0394 BCAA (kJ/g)")
pred_titles_pt2[[2]]<-"Control Diet BCAA (kJ/g)" 
pred_titles_pt2[[3]]<-"ma"
pred_titles_pt2[[4]]<-"Control Diet P to C Ratio" 
pred_titles_pt2[[5]]<-"Control Diet Protein (kJ/g)" 

pred_titles_pt3<-list()
pred_titles_pt3[[1]]<-"Control Diet BCAA (kJ/g)"
pred_titles_pt3[[2]]<-c("Control Diet P to C Ratio", "\U0394 BCAA (kJ/g)") 
pred_titles_pt3[[3]]<-c("\U0394 BCAA to non-BCAA ratio (kJ/g)")
pred_titles_pt3[[4]]<-c("\U0394 BCAA to non-BCAA ratio (kJ/g)") 

pred_titles<-list(pred_titles_pt1, pred_titles_pt2, pred_titles_pt3)

########################### Orchard Plots for MS

MA_plot_list<-list()

# For traits in each 'part'
for(p in 1:length(traits)){
	
	# Pull out the pth set of results and pth set of data	
	models_p<-MA_pool_list[[p]]
	data_p<-data_list[[p]]
	TSV_p<-TSV_list[[p]]
	traits_p<-traits[[p]]
	
	for(t in 1:length(models_p)){
		
		# Get the tth model
		model_t<-models_p[[t]][[1]]
		
		# Get the tth set of MA results
		ma_res_t<-data.frame(name=traits_p[t], estimate=model_t[[1]][,1])
		ma_res_t$lowerCL<-ma_res_t$estimate - 1.96 * sqrt(model_t[[1]][,2])
		ma_res_t$upperCL<-ma_res_t$estimate + 1.96 * sqrt(model_t[[1]][,2])
		ma_res_t$lowerPR<-ma_res_t$estimate - 1.96 * sqrt(sum(model_t[[2]][,c(1,2)]) + TSV_p[[t]][[1]])
		ma_res_t$upperPR<-ma_res_t$estimate + 1.96 * sqrt(sum(model_t[[2]][,c(1,2)]) + TSV_p[[t]][[1]])
		
		# Get the tth set of set
		data_orchard_t<-data_p[[t]][[1]]
		
		# Bind together
		if(t == 1){
			ma_res<-ma_res_t
			data_orchard<-data_orchard_t
		}else{
			ma_res<-rbind(ma_res, ma_res_t)
			data_orchard<-rbind(data_orchard, data_orchard_t)
		}
		
	}
	
	# Create the orchard object
	orchard_object<-list(mod_table=ma_res, data=data_orchard)
	
	# Plot it out
	MA_plot_list[[p]]<-my_orchard(object=orchard_object, xlab="lnRR", mod="Trait", mod_labels=trait_titles[[p]], alpha=0.25)
	
}



########################### Orchard Plots for Fast vs fed

MR_fast_fed<-list()

# For traits in each 'part'
for(p in 1:length(traits[-3])){
	
	# Pull out the pth set of results and pth set of data	
	models_p<-MR_fast_list[[p]]
	data_p<-data_list[[p]]
	TSV_p<-TSV_list[[p]]
	traits_p<-traits[[p]]
	
	for(t in 1:length(models_p)){
		
		# Get the tth model
		model_t<-models_p[[t]]
		
		# Get the tth set of MA results
		ma_res_t<-data.frame(name=paste0(traits_p[t], "_", row.names(model_t[[1]])), estimate=model_t[[1]][,1])
		ma_res_t$lowerCL<-ma_res_t$estimate - 1.96 * sqrt(diag(model_t[[1]][,-1]))
		ma_res_t$upperCL<-ma_res_t$estimate + 1.96 * sqrt(diag(model_t[[1]][,-1]))
		ma_res_t$lowerPR<-ma_res_t$estimate - 1.96 * sqrt(sum(model_t[[2]][,c(1,2)]) + TSV_p[[t]][[1]])
		ma_res_t$upperPR<-ma_res_t$estimate + 1.96 * sqrt(sum(model_t[[2]][,c(1,2)]) + TSV_p[[t]][[1]])
		
		# Get the tth set of set
		data_orchard_t<-data_p[[t]][[1]]
		data_orchard_t$Pred<-paste0(data_orchard_t$Trait, "_", "Time_fasted", data_orchard_t$Time_fasted)
		
		# Bind together
		if(t == 1){
			ma_res<-ma_res_t
			data_orchard<-data_orchard_t
		}else{
			ma_res<-rbind(ma_res, ma_res_t)
			data_orchard<-rbind(data_orchard, data_orchard_t)
		}
		
	}
	
	# Create the orchard object
	orchard_object<-list(mod_table=ma_res, data=data_orchard)
	
	# Create new nice titles
	titles<-NA
	for(g in 1:length(trait_titles[[p]])){
		titles<-c(titles, paste0(trait_titles[[p]][g], "\nFasted"), paste0(trait_titles[[p]][g], "\nFed"))
	}
	titles<-titles[-1]
	
	# Plot it out
	MR_fast_fed[[p]]<-my_orchard(object=orchard_object, xlab="lnRR", mod="Pred", mod_labels=titles, alpha=0.25, rep_cb=2)
	
}


################################## Bubble plots and surfaces for nutritional predictors
 
MR_plot_list<-list()
surface_significance_list<-list()
signif_count<-1

for(p in 1:length(traits)){
	
	# Plot list for part p
	plot_list_p<-list()
	
	# Fot each trait in the part
	for(t in 1:length(traits[[p]])){
		
		print(traits[[p]][t])
		
		# Pull out the model, model type and dataset for pt
		model_t<-aic_models[[p]][[t]][[1]]
		type_t<-aic_types[[p]][[t]]
		data_t<-data_list[[p]][[t]][[1]]
		trait_t<-traits[[p]][[t]]
		pred_t<-pred_titles[[p]][[t]]
		
		# If we have a 1d linear do a bubble plot
		if(type_t == "lin.1d"){
			moderator<-row.names(model_t)[2]
			plot_list_p[[t]]<-my_bubble(model=model_t, data=data_t, moderator=moderator, xlab=pred_t, ylab=paste0(trait_titles[[p]][[t]], " (lnRR)"), alpha=0.25)
		}
		
		# If we have a 1d spline also do the bubble plot, but with the spline
		if(type_t == "spline.1d"){
			moderator<-row.names(model_t)[2]
			moderator<-strsplit(moderator, "(", fixed=T)[[1]][2]
			moderator<-strsplit(moderator, ",", fixed=T)[[1]][1]
			plot_list_p[[t]]<-my_bubble(model=model_t, data=data_t, moderator=moderator, xlab=pred_t, ylab=paste0(trait_titles[[p]][[t]], " (lnRR)"), alpha=0.25, spline=T, df=df)
		}
		
		# If we have a 2d linear do a surface
		if(type_t == "lin.2d"){
			moderator1<-row.names(model_t)[2]
			moderator2<-row.names(model_t)[3]
			surfaces<-my_surface(model=model_t, data=data_t, mods=c(moderator1, moderator2), xlab=pred_t[1], ylab=pred_t[2], title=paste0(trait_titles[[p]][[t]], " (lnRR)"), lab_size=3)
			plot_list_p[[t]]<-surfaces[[1]]
			surface_significance_list[[signif_count]]<-surfaces[[2]]
			signif_count<-signif_count + 1			
		}	
		
		# We do not have 2d splines, but if we did, we could add them in here
		if(type_t == "spline.2d"){
			moderator1<-row.names(model_t)[2]
			moderator1<-strsplit(moderator1, "(", fixed=T)[[1]][2]
			moderator1<-strsplit(moderator1, ",", fixed=T)[[1]][1]
			moderator2<-row.names(model_t)[2]
			moderator2<-strsplit(moderator2, ":", fixed=T)[[1]][2]
			moderator2<-strsplit(moderator2, "(", fixed=T)[[1]][2]
			moderator2<-strsplit(moderator2, ",", fixed=T)[[1]][1]
			surfaces<-my_surface(model=model_t, data=data_t, mods=c(moderator1, moderator2), xlab=pred_t[1], ylab=pred_t[2], title=paste0(trait_titles[[p]][[t]], " (lnRR)"), lab_size=3, spline=T, df=df)
			plot_list_p[[t]]<-surfaces[[1]]
			surface_significance_list[[signif_count]]<-surfaces[[2]]
			signif_count<-signif_count + 1
		}
		
		# If its a MA just add NA to the plot list we wont be plotting anything
		if(type_t == "ma"){
			plot_list_p[[t]]<-NA
		}
		
	}
	
	# Save the pth set of plots
	MR_plot_list[[p]]<-plot_list_p

}

###################### Lets start laying some of these out

tag_size<-25

# Figure 1 - circulating levels of BCAAs
CairoPDF("Figure_1.pdf", height=5*2, width=5*4)

# Grab the right plots and add a few extras for the arrangement
A<-MA_plot_list[[1]] + ggtitle("A.") + theme(plot.title=element_text(size=25, face="bold")) + theme(plot.margin=margin(t=1.5,r=20,b=1.5,l=1.5))
B<-MR_fast_fed[[1]] + ggtitle("B.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))
C<-MR_plot_list[[1]][[1]] + ggtitle("C.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))
D<-MR_plot_list[[1]][[2]] + ggtitle("D.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))
E<-MR_plot_list[[1]][[3]] + ggtitle("E.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))
Fa<-MR_plot_list[[1]][[4]] + ggtitle("F.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))

grid.arrange(A, B, C, D, E, Fa, layout_matrix=rbind(c(1,1,1,2,2,2,3,3,3,5,5,5),
												    c(1,1,1,2,2,2,3,3,3,5,5,5),
												    c(1,1,1,2,2,2,3,3,3,5,5,5),
												    c(1,1,1,2,2,2,4,4,4,6,6,6),
												    c(1,1,1,2,2,2,4,4,4,6,6,6),
												    c(1,1,1,2,2,2,4,4,4,6,6,6)))



dev.off()

# Figure 2 Glucose Homeostasis
CairoPDF("Figure_2.pdf", height=5*2, width=5*4)

# Grab the right plots and add a few extras for the arrangement
A<-MA_plot_list[[2]] + ggtitle("A.") + theme(plot.title=element_text(size=25, face="bold")) + theme(plot.margin=margin(t=1.5,r=20,b=1.5,l=1.5))
B<-MR_fast_fed[[2]] + ggtitle("B.") + theme(plot.title=element_text(size=25, face="bold")) + theme(plot.margin=margin(t=1.5,r=20,b=1.5,l=1.5))
C<-MR_plot_list[[2]][[1]] + ggtitle("C.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))
D<-MR_plot_list[[2]][[2]] + ggtitle("D.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))
E<-MR_plot_list[[2]][[4]] + ggtitle("E.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))
Fa<-MR_plot_list[[2]][[5]] + ggtitle("F.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))


grid.arrange(A, B, C, D, E, Fa, layout_matrix=rbind(c(1,1,1,2,2,2,3,3,3,5,5,5),
												    c(1,1,1,2,2,2,3,3,3,5,5,5),
												    c(1,1,1,2,2,2,3,3,3,5,5,5),
												    c(1,1,1,2,2,2,4,4,4,6,6,6),
												    c(1,1,1,2,2,2,4,4,4,6,6,6),
												    c(1,1,1,2,2,2,4,4,4,6,6,6)))

dev.off()

# Figure 3 Body Composition and Food Intake
CairoPDF("Figure_3.pdf", height=5*2, width=5*3)

# Grab the right plots and add a few extras for the arrangement
A<-MA_plot_list[[3]] + ggtitle("A.") + theme(plot.title=element_text(size=25, face="bold")) + theme(plot.margin=margin(t=1.5,r=20,b=1.5,l=1.5))
B<-MR_plot_list[[3]][[1]] + ggtitle("B.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))
C<-MR_plot_list[[3]][[2]] + ggtitle("C.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))
D<-MR_plot_list[[3]][[3]] + ggtitle("D.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))
E<-MR_plot_list[[3]][[4]] + ggtitle("E.") + theme(plot.title=element_text(size=tag_size, face="bold")) + theme(plot.margin=margin(t=1.5,r=10,b=1.5,l=1.5))

grid.arrange(A, B, C, D, E, layout_matrix=rbind(c(1,1,1,2,2,2,4,4,4),
												c(1,1,1,2,2,2,4,4,4),
												c(1,1,1,2,2,2,4,4,4),
												c(1,1,1,3,3,3,5,5,5),
												c(1,1,1,3,3,3,5,5,5),
												c(1,1,1,3,3,3,5,5,5)))

dev.off()

# Figure significnace of surfaces - to be embedded usuing external graphic software
CairoPDF("surface_signif.pdf", height=5, width=5)
for(i in 1:length(surface_significance_list)){
	print(surface_significance_list[[i]])
}

dev.off()

