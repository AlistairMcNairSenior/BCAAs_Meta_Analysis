
# Function to get typical sampling variance by AM Senior @ the university of Sydney on the 26.sept.2016

my_TMV<-function(Sampling.Variance){

	MEV<-((sum(1/Sampling.Variance)) * (length(Sampling.Variance)-1)) / (((sum(1/Sampling.Variance))^2) - (sum((1/Sampling.Variance)^2)))

return(MEV)

}

# Function for aggregating m models from rma.mv from m imputed datasets - the models must all have the same number of dimensions - takes a list, of length m and returns the aggregate estimates and VCOV

my_agg<-function(model_list){
	
	# Check we have rma.mv models
	screen<-lapply(model_list, class) == "data.frame"
	tag<-which(screen == T)
	if(length(tag) < length(model_list)){
		
		# Drop any failues from the 
		if(length(tag) > 0){
			for(i in 1:length(tag)){
				location<-tag[i] - (i-1)
				model_list[[location]]<-NULL
			}
		}
		# Coefficients in the model
		coef<-row.names(model_list[[1]]$b)
		n_coef<-length(coef)
		n_ranef<-length(model_list[[1]]$sigma2)
		use_tau<-F
		if(n_ranef == 0){
			use_tau<-T
			n_ranef<-1
		}
		
		# Number of imputes
		m<-length(model_list)
		
		# To hold the pooled values
		# point estimates
		Q_m<-t(as.matrix(array(0, c(1, n_coef))))
		# Pooled within VCV
		U_m<-as.matrix(array(0, c(n_coef, n_coef)))
		# Between VCV
		B_m<-as.matrix(array(0, c(n_coef, n_coef)))
		# variance estimates
		sigma2<-as.matrix(array(0, c(1, n_ranef)))
		# Q_between for moderators and the p value therefore
		QB<-as.matrix(array(0, c(1, 2)))
		aic<-0
		
		# Do the average for each
		for(i in 1:m){
			Q_m<-Q_m + model_list[[i]]$b / m
			U_m<-U_m + vcov(model_list[[i]]) / m
			if(use_tau == T){
				sigma2<-sigma2 + model_list[[i]]$tau2/ m
			}else{
				sigma2<-sigma2 + model_list[[i]]$sigma2/ m
			}
			QB<-QB + c(model_list[[i]]$QM, model_list[[i]]$QMp) / m
			aic<-aic + AIC(model_list[[i]]) / m
		}
		
		# Get the between impute variance
		for(i in 1:m){
			B_m<-B_m + ((model_list[[i]]$b - Q_m) %*% t(model_list[[i]]$b - Q_m)) / (m-1)
		}
		
		# The total Variance
		T_m<-U_m + (1 + m^-1) * B_m
		
		# Organise in to a martrix where column one is the betas, and subsequent columns are the VCV
		out<-list()
		out[[1]]<-cbind(Q_m, T_m)
		row.names(out[[1]])<-coef
		colnames(out[[1]])<-c("b", paste0("v_", coef))
		out[[2]]<-t(as.matrix(c(sigma2, QB)))
		colnames(out[[2]])<-c(paste0("sigma2_", seq(1, n_ranef, 1)), c("QM", "QMp"))
		
		# Check we have an rma.uni.trimfill model and add a bit more info
		if(class(model_list[[1]])[1] == "rma.uni.trimfill"){
			
			k_missing<-0
			n_right<-0
			for(i in 1:m){
				k_missing<-k_missing + model_list[[i]]$k0 / m
				n_right<-n_right + as.numeric((model_list[[i]]$side == "right")) / m
			}
			
			# Add in the info on the missing studies
			out[[3]]<-c(k_missing, n_right)
			names(out[[3]])<-c("Missing Effects", "prop_right")
		}else{
			out[[3]]<-aic
			names(out[[3]])<-"AIC"
		}
		out[[4]]<-paste0("N models = ", m)
	}else{
		
		# If there is no model return a bunch of NAs
		out<-list()
		out[[1]]<-as.matrix(array(NA, c(1, 2)))
		row.names(out[[1]])<-"NoModel"
		colnames(out[[1]])<-c("b_NoModel", paste0("v_NoModel"))
		out[[2]]<-NA
		out[[3]]<-NA
		
	}
	
	# Return the output	
	return(out)
}


# Function to get pairwise data
my_pairwise<-function(data, col_unit="Experimental_Unit"){
		
		exps<-unique(data[, col_unit])
		
		for(i in 1:length(exps)){
		
			# Pull out the ith experrimental unit
			exp_i<-data[which(data[, col_unit] == exps[i]),]
			
			# Ensure it is sorted by total dietary BCAAs (should already be)
			exp_i<-exp_i[order(exp_i$Diet_BCAA_kJ_g),]
			
			# Get all the pairwise combinations of diets within experiment
			combinations<-t(combn(exp_i$Group, m=2))
			
			# Pull out the correct data and bind together pairwise
			control<-exp_i[match(combinations[,1], exp_i$Group),]
			
			targets<-c("Group", "n", "Mean", "SD")
			targets<-c(grep("Diet", names(exp_i)), match(targets, names(exp_i)))
			
			group<-exp_i[match(combinations[,2], exp_i$Group), targets]
			names(group)<-paste0("exp_", names(group))
			pw<-as.data.frame(cbind(control, group))
			delta_BCAA<-pw$exp_Diet_BCAA_kJ_g - pw$Diet_BCAA_kJ_g
			
			# Double check we have no comparisons among identical BCAA levels
			drop<-which((pw$delta_BCAA) <= 0)
			if(length(drop >= 1)){
				pw<-pw[-drop,]
			}
			
			# Add to the pairwise dataset
			if(i == 1){
				pw_data<-pw
			}else{
				pw_data<-rbind(pw_data, pw)	
			}
			
		}
		
		pw_data$ES_ID<-paste0(pw_data$Group, "x", pw_data$exp_Group)
		return(pw_data)
}


# A function created to find the outer perimeter over which the surface should be fitted
findConvex<-function(x,y,rgnames,res=101, x.new=NA, y.new=NA){
	hull<-cbind(x,y)[chull(cbind(x,y)),]
	px<-pretty(x)
	py<-pretty(y)
	if(is.na(x.new[1])==TRUE){
		x.new<-seq(min(px, na.rm=T),max(px, na.rm=T),len=res)
	}
	if(is.na(y.new[1])==TRUE){
		y.new<-seq(min(py, na.rm=T),max(py, na.rm=T),len=res)
	}
	ingrid<-as.data.frame(expand.grid(x.new,y.new))                                                              
	Fgrid<-ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
	names(Fgrid)<-rgnames
	return(Fgrid)
}


my_orchard<-function(object, mod, xlab, mod_labels=NA, mod_title="", alpha = 0.5, angle = 90, cb = TRUE, rep_cb = 1) 
{
	
	require(plyr)
	
	data <- object$data
	
	data$moderator<-data[,mod]
	
	data$scale <- (1/sqrt(data[, "vi"]))
	
	legend <- "Precision (1/SE)"
	      
	label <- xlab
	
	if(is.na(mod_labels[1]) == T){
		mod_labels<-object$mod_table$name
	}
	
	# Add in the Ks    
	ks<-ddply(data, .(moderator), summarise, length(yi))    
	object$mod_table$K <- ks$..1[match(object$mod_table$name, ks$moderator)]
	
	# Make sure we have the factors ordered as they are in the table
	object$mod_table$name<-factor(object$mod_table$name, levels=object$mod_table$name)
	data$moderator<-factor(data$moderator, levels=object$mod_table$name)
	
	cbpl <- c("#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#56B4E9", "#999999")
	cbpl <- matrix(apply(matrix(cbpl, ncol=1), 1, rep, rep_cb), nrow=1)
	
	plot <- ggplot2::ggplot(data = object$mod_table, aes(x = estimate, y = name)) + 
		ggbeeswarm::geom_quasirandom(data = data, aes(x = yi, y = moderator, size = scale, colour = moderator), groupOnX = FALSE, alpha = alpha) + 
		ggplot2::geom_errorbarh(aes(xmin = object$mod_table$lowerPR, xmax = object$mod_table$upperPR), height = 0, show.legend = FALSE, size = 0.5, alpha = 0.6) + 
		ggplot2::geom_errorbarh(aes(xmin = object$mod_table$lowerCL, xmax = object$mod_table$upperCL), height = 0, show.legend = FALSE, size = 1.2) + 
		ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
		ggplot2::geom_point(aes(fill = object$mod_table$name), size = 3, shape = 21) + 
		ggplot2::annotate("text", x = (max(data$yi) + (max(data$yi) * 0.1)), y = (seq(1, dim(object$mod_table)[1], 1) + 0.3), label = paste("italic(k)==", object$mod_table$K), parse = TRUE, hjust = "right", size = 5) + 
		ggplot2::theme_bw() + 
		ggplot2::guides(fill = "none", colour = "none", size=guide_legend(title.position="top", title.hjust=0.5)) + 
		ggplot2::theme(legend.position = c(1, 0), legend.justification = c(1, 0)) + 
		ggplot2::theme(legend.title = element_text(size = 10)) +
		ggplot2::theme(legend.direction = "horizontal") + 
		ggplot2::theme(legend.background = element_blank()) +         
		ggplot2::labs(x = label, y = "", size = legend) +
		ggplot2::scale_y_discrete(labels=mod_labels, name=mod_title) +
		ggplot2::theme(axis.text.y = element_text(size = 15, colour = "black", hjust = 0.5, angle = angle), axis.text.x = element_text(size=15), axis.title.x = element_text(size=15))
	
	if (cb == TRUE) {
		plot <- plot + scale_fill_manual(values = cbpl) + scale_colour_manual(values = cbpl)
	}
	    		
	return(plot)
}

# Function to make bubble plots
my_bubble<-function(model, data, xlab=NA, ylab, moderator, alpha=0.5, est_col="red", est_width=1, spline=F, df=3){
	
	# Load the libraries
	require(ggplot2)
	require(splines)
	
	# Set the xlabel as the moderator if its missing
	if(is.na(xlab) == T){
		xlab<-moderator
	}
		
	# Add in the right moderator column and the scale
	data$moderator<-data[,moderator]
	data$scale<-1/sqrt(data$vi)
	
	# Get the X values to predict over
	pred<-seq(min(data[,moderator]), max(data[,moderator]), length=100)
	
	# Get the spline if we are doing that
	if(spline == T){
		sp1<-bs(data[,moderator], df=df)
		X<-model.matrix(~ predict(sp1, pred))
	}else{
		X<-model.matrix(~ pred)
	}
	
	
	# Get the estimate as Y and the SE
	Y<-X %*% model[,1]
	SE<-sqrt(diag(X %*% model[,-1] %*% t(X)))
	
	# Aggregate the data	
	plot_data<-data.frame(pred=pred, Y=Y, LCL=Y-1.96*SE, UCL=Y+1.96*SE)
	
	# Make the plot
	plot<-ggplot2::ggplot(data=plot_data, aes(x=pred, y=Y)) +
		ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black") +
		ggplot2::geom_point(data=data, aes(x=moderator, y=yi, size=scale), alpha=alpha, col="#0072B2") +
		ggplot2::geom_ribbon(aes(x=pred, ymin=LCL, ymax=UCL), fill=est_col, alpha=alpha) +
		ggplot2::geom_line(aes(x=pred, y=Y), color=est_col, size=est_width) +
		ggplot2::theme_bw() + 
		ggplot2::guides(fill = "none", colour = "none", size=guide_legend(title.position="top", title.hjust=0.5, override.aes=list(color="black"))) + 
		ggplot2::theme(legend.position = c(1, 0), legend.justification = c(1, 0)) + 
		ggplot2::theme(legend.title = element_text(size = 10)) +
		ggplot2::theme(legend.direction = "horizontal") +
		ggplot2::labs(x = xlab, y = ylab, size="Precision (1/SE)") + 
		ggplot2::theme(legend.background = element_blank()) +
		ggplot2::theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size=15), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))       
	
	return(plot)
}


# Function to make the surfaces
my_surface<-function(model, data, mods, res=301, title="", xlab="", ylab="", add_effects=F, max_size=5, nlevels=5, contour_at=NA, x_lim=NA, y_lim=NA, lab_size=5, spline=F, df=3){
	
	# Add in predictor columns
	data$pred1<-data[,mods[1]]
	data$pred2<-data[,mods[2]]
		
	# Pull the coefficients out of the model
	coef<-model[,1]
	
	# Pretty range over which to fit
	new_x<-seq(min(pretty(data$pred1)), max(pretty(data$pred1)), length=res)
	new_y<-seq(min(pretty(data$pred2)), max(pretty(data$pred2)), length=res)
	
	# Values to predict over
	pred_vals<-findConvex(data$pred1, data$pred2, c("pred1", "pred2"), res=res, x.new=new_x, y.new=new_y)
	pred_vals<-pred_vals[which(apply((is.na(pred_vals) == F), 1, prod) == 1),]
	
	# Get the model matrix
	# If we are doing splines get those
	if(spline == T){
		sp1<-bs(data$pred1, df=df)
		sp2<-bs(data$pred2, df=df)
		X<-model.matrix(~ predict(sp1, pred_vals[,1]) : predict(sp2, pred_vals[,2]))
	}else{
		X<-model.matrix(~ pred_vals[,1] * pred_vals[,2])
	}
	
	# Do the prediction
	Y<-X %*% coef
	SE<-sqrt(diag(X %*% model[,-1] %*% t(X)))
	
	# Define significance
	CI<-abs(Y) - 1.96 * SE
	signif<-rep(0, length(CI))
	signif[which(CI > 0)]<-1
	
	# Find the minimum and maximum values of the surface
	mx<-max(abs(range(Y, na.rm=TRUE)))
	mn<-mx * -1
	locs<-(range(Y, na.rm=TRUE) - mn) / (mx-mn) * 256
	
	# Get the colors to use from the pallette specified above
	rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
	map<-rgb.palette(256)
	
	# Arrange the data to plot
	plot_data<-cbind(pred_vals, Y, signif)
	points_data<-data.frame(pred1=data$pred1, pred2=data$pred2)
	points_data$col_code<-rep("red", length(data$yi))
	points_data$col_code[which(data$yi < 0)]<-"blue"
	points_data$col_code[which((abs(data$yi) - 1.96*sqrt(data$vi)) < 0)]<-"grey"
	points_data$size_code<-1/sqrt(data$vi)
	points_data$size_code<-points_data$size_code / max(points_data$size_code) * max_size
	if(add_effects == F){
		points_data$size_code<--1
	}
	
	# Set the axes to default to the estimated rangess
	if(is.na(x_lim[1]) == T){
		x_lim=range(new_x)
	}
	if(is.na(y_lim[1]) == T){
		y_lim=range(new_y)
	}
	
	# Set the contour
	if(is.na(contour_at) == T){
		contour_at<-signif(seq(0, mx, length=nlevels/2)[2], 1)
	}
	
	# Make the surface in ggplot2	
	plot<-ggplot2::ggplot(data=plot_data, aes(x=pred1, y=pred2)) +
		ggplot2::geom_raster(aes(fill=Y), show.legend=F, interpolate=F, na.rm=T) +
		ggplot2::scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
		ggplot2::geom_contour(data=plot_data, aes(x=pred1, y=pred2, z=Y), na.rm=T, color="black", binwidth=contour_at) +
		metR::geom_label_contour(data=plot_data, aes(x=pred1, y=pred2, z=Y), size=lab_size, binwidth=contour_at) +
		#ggplot2::geom_abline(slope=4, intercept=0, linetype = 2, colour = "purple") +
		ggplot2::labs(x = xlab, y = ylab) +
		ggplot2::xlim(x_lim) +
		ggplot2::ylim(range(y_lim)) +
		ggplot2::theme_bw() +
		ggplot2::geom_point(data=points_data, aes(x=pred1, y=pred2), na.rm=T, color=points_data$col_code, size=points_data$size_code, shape=1) +
		ggplot2::theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size=15), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15)) +
		ggplot2::annotate("text", x = min(x_lim), y = max(y_lim), label = title, hjust = 0, vjust = 1, size = 5)
	
	# Make a plot of the significance
	plot2<-ggplot2::ggplot(data=plot_data, aes(x=pred1, y=pred2)) +
			 ggplot2::geom_raster(aes(fill=signif), show.legend=F, interpolate=F, na.rm=T) +
			 ggplot2::labs(x = NULL, y = NULL) +
			 ggplot2::theme_bw() +
			 ggplot2::theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(),
			 				axis.ticks.y=element_blank(), axis.text.y=element_blank()) +				
			 ggplot2::scale_fill_gradientn(colors=c("grey", "purple")) +
			 ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", size = 1.5) +
			 ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "black", size = 1.5)

	
	# Return the plots
	plots_list<-list(plot, plot2)
	return(plots_list)
	
}

# Do all the analyses for the traits specified
analyse_traits<-function(traits, nutri_formulas=list(), nutri_params=list(), other_formulas=list(), other_params=list(), n_imp=1, wd, core=1){
	
	# Load the tcltk package for the progress bar
	library(mice)
	library(metafor)
	library(plyr)
	library(corpcor)
	library(ape)
	library(splines)
	library(tcltk)
	
	pb<-tkProgressBar(min = 0, max = (length(traits) * (length(nutri_formulas) + length(other_formulas)) + length(traits)) * n_imp, title = paste0("Progress on Core ", core))
	progress<-0
	
	for(t in 1:length(traits)){	
		
		# Pull out the right trait
		trait_t<-traits[t]
		
		# Read in data for the trait
		setwd(paste0(wd, "/data"))
		data<-read.csv(paste0("final_", trait_t, ".csv"))
		
		# Double check that we do not have missing nutritional data and drop any if we do
		test<-apply(data[, c("Diet_P", "Diet_C", "Diet_F", "Diet_BCAA", "Diet_kJ_g")], 2, is.na)
		tag<-which(apply(test, 1, sum) > 0)
		if(length(tag) > 0){
			data<-data[-tag,]
		}
		
		# Add some columns for the absolute amount of energy coming from nutrients and also ratios
		data$Diet_P_kJ_g<-data$Diet_P/100 * data$Diet_kJ_g
		data$Diet_C_kJ_g<-data$Diet_C/100 * data$Diet_kJ_g
		data$Diet_F_kJ_g<-data$Diet_F/100 * data$Diet_kJ_g
		data$Diet_BCAA_kJ_g<-data$Diet_BCAA/100 * data$Diet_kJ_g
		data$Diet_nonBCAA_kJ_g<-data$Diet_P_kJ_g - data$Diet_BCAA_kJ_g
		data$Diet_Iso_kJ_g<-data$Diet_Iso/100 * data$Diet_kJ_g
		data$Diet_Leu_kJ_g<-data$Diet_Leu/100 * data$Diet_kJ_g
		data$Diet_Val_kJ_g<-data$Diet_Val/100 * data$Diet_kJ_g
		data$Diet_PC_kJ_g<-data$Diet_P_kJ_g / data$Diet_C_kJ_g
		data$Diet_PF_kJ_g<-data$Diet_P_kJ_g / data$Diet_F_kJ_g
		data$Diet_BCAA.Non_BCAA_kJ_g<-data$Diet_BCAA_kJ_g / data$Diet_nonBCAA_kJ_g
		# Convert Time_Fatsed to fasted vs fed
		Fast_Fed<-rep("Fasted", dim(data)[1])
		tag<-which(data$Time_fasted == 0)
		if(length(tag) > 0){
			Fast_Fed[tag]<-"Fed"
		}
		rm(tag)
		tag<-which(is.na(data$Time_fasted) == T)
		if(length(tag) > 0){
			Fast_Fed[tag]<-NA
		}	
		data$Time_fasted<-Fast_Fed
				
		# Do the imputation for the missing SDs, if necessary
		missing<-which(is.na(data$SD) == T)
		if(length(missing) > 0){
			log_vals<-data.frame(ln_mu=log(abs(data$Mean)), ln_SD=log(data$SD))
			my_imp<-mice(log_vals, m=n_imp, print=FALSE)$imp$ln_SD
			my_imp<-exp(my_imp)
		}
				
		# Results objects for trait t
		trait_t_MAs<-list()
		trait_t_TSVs<-data.frame(rep=c(1:n_imp), TSV=0)
		trait_t_MRs<-as.list(rep(0, length(nutri_formulas)))
		names(trait_t_MRs)<-nutri_formulas
		trait_t_MRs<-lapply(trait_t_MRs, as.list)
		trait_t_MRs_2<-as.list(rep(0, length(other_formulas)))
		names(trait_t_MRs_2)<-other_formulas
		trait_t_MRs_2<-lapply(trait_t_MRs_2, as.list)
		trait_t_ERs<-list()
		trait_t_TFs<-list()
		
		# Record if we are bending the VCV matrix
		bend<-0
		
		# Doing the analysis n_imp times, with any imputed values
		for(m in 1:n_imp){
			
			# Add in the imputed values
			data_m<-data
			if(length(missing) > 0){
				data_m$SD[as.numeric(row.names(my_imp))]<-my_imp[,m]
			}
			
			# Format pairwise comparisons within the experiment
			pw_data<-my_pairwise(data=data_m)
			
			# Calculate the effects sizes
			pw_data<-escalc(n1i=exp_n, n2i=n, m1i=exp_Mean, m2i=Mean, sd1i=exp_SD, sd2i=SD, measure="ROM", data=pw_data)
			
			# Calculate the difference in dietary levels for meta-regression
			pw_data$diff_P<-pw_data$exp_Diet_P_kJ_g - pw_data$Diet_P_kJ_g
			pw_data$diff_BCAA<-pw_data$exp_Diet_BCAA_kJ_g - pw_data$Diet_BCAA_kJ_g
			pw_data$diff_nonBCAA<-pw_data$exp_Diet_nonBCAA_kJ_g - pw_data$Diet_nonBCAA_kJ_g
			pw_data$diff_Iso<-pw_data$exp_Diet_Iso_kJ_g - pw_data$Diet_Iso_kJ_g
			pw_data$diff_Leu<-pw_data$exp_Diet_Leu_kJ_g - pw_data$Diet_Leu_kJ_g
			pw_data$diff_Val<-pw_data$exp_Diet_Val_kJ_g - pw_data$Diet_Val_kJ_g
			pw_data$diff_BCAA.Non_BCAA_kJ_g<-pw_data$exp_Diet_BCAA.Non_BCAA_kJ_g - pw_data$Diet_BCAA.Non_BCAA_kJ_g
			
			# Calculate the TSV
			trait_t_TSVs$TSV[m]<-my_TMV(pw_data$vi)

			# Calculate the VCV matrix - a real ball ache here
								
			# Create the VCV for effect sizes with shared data
			VCV<-as.matrix(array(0, c(dim(pw_data)[1],dim(pw_data)[1])))
			rownames(VCV)<-pw_data$ES_ID
			colnames(VCV)<-pw_data$ES_ID
			
			# flag any control groups with more than one entry, and any treatment groups with more than one entry, or that also appear in the control column
			low.covs<-ddply(pw_data, .(Group), here(summarise), flag=length(n) > 1)
			high.covs<-ddply(pw_data, .(exp_Group), here(summarise), flag1=(length(n) > 1), flag2=(sum(pw_data$Group == exp_Group[1]) >= 1))
			
			# Go through the low flags and add in the covariances	
			for(b in 1:length(low.covs$Group)){
				# If it has been flagged
				if(low.covs$flag[b] == T){
					# Look for the ID in the low IDs and get the combinations of rows / columns
					coord<-which(pw_data$Group == low.covs$Group[b])
					coord<-combn(coord, m=2)
					for(v in 1:dim(coord)[2]){
						# Go through and add in
						cov_v<-pw_data$SD[coord[1,v]]^2 / (pw_data$n[coord[1,v]] * pw_data$Mean[coord[1,v]]^2)
						VCV[coord[1,v], coord[2,v]]<-cov_v
						VCV[coord[2,v], coord[1,v]]<-cov_v
					}	
				}
			}
			
			# Repeat for the high flags
			for(b in 1:length(high.covs$exp_Group)){
				if(high.covs$flag1[b] == T | high.covs$flag2[b] == T){
					coord<-which(pw_data$exp_Group == high.covs$exp_Group[b] | pw_data$Group == high.covs$exp_Group[b])
					coord<-combn(coord, m=2)
					for(v in 1:dim(coord)[2]){
						cov_v<-pw_data$exp_SD[coord[1,v]]^2 / (pw_data$exp_n[coord[1,v]] * pw_data$exp_Mean[coord[1,v]]^2)
						VCV[coord[1,v], coord[2,v]]<-cov_v
						VCV[coord[2,v], coord[1,v]]<-cov_v
					}	
				}
			}
			
			# Add in the diagonal elements
			diag(VCV)<-pw_data$vi
			
			# Bend if necessary
			if(is.positive.definite(VCV) == F){
				VCV<-make.positive.definite(VCV)
				bend<-bend + 1
			}
			
			# Fit the MA
			RMA<-rma.mv(yi=yi, V=VCV, random=list(~1|Experimental_Unit, ~1|ES_ID), data=pw_data)
			
			# Save the results for the mth imputation
			trait_t_MAs[[m]]<-RMA
			
			# PB tests on the MA
			# Make a residual that removes the experimental estimate
			prediction<-RMA$b[1] + ranef(RMA)$Experimental_Unit[match(pw_data$Experimental_Unit, row.names(ranef(RMA)$Experimental_Unit)),1] 
			my_res<-pw_data$yi - prediction
			
			# Do some PB tests
			pb_mod<-rma(yi=resid(RMA), vi=pw_data$vi)
			trait_t_ERs[[m]]<-regtest(pb_mod)$fit
			trait_t_TFs[[m]]<-trimfill(pb_mod)
			
			# Update the progress
			progress<-progress + 1
			setTkProgressBar(pb, progress)
			
			# Run the nutritional formuals for which there are sufficient data						
			for(r in 1:length(nutri_formulas)){
				
				# Test that there are sufficient data
				test_n<-(nutri_params[[r]] * 10) <= dim(pw_data)[1]
				
				# run the meta-regression if there are enough data - otherwise just record NAs
				if(test_n == T){
					# for a numeric predictor fit a regular regression
					form<-paste0("~", nutri_formulas[[r]])
					model<-try(rma.mv(yi=yi, V=VCV, random=list(~1|Experimental_Unit, ~1|ES_ID), data=pw_data, mods=as.formula(form)), silent=T)				
				# Otherwise set the model as a dataframe
				}else{
					model<-data.frame(1)
				}
				# Also if the model fails set as a dataframe
				if(class(model)[1] == "try-error"){
					model<-data.frame(1)
				}

				# Save the reuslts for the rth regression in the mth imputation	
				trait_t_MRs[[r]][[m]]<-model
				
				# Clean up a bit
				rm(model)
				
				# Update the progress
				progress<-progress + 1
				setTkProgressBar(pb, progress)
			}
			
			# Run the other formulas for which there are sufficient data						
			for(r in 1:length(other_formulas)){
				
				# Check and remove missing predictors
				drop<-which(is.na(pw_data[, other_formulas[[r]]]) == T)
				pw_data_r<-pw_data
				VCV_r<-VCV
				if(length(drop) > 0){
					pw_data_r<-pw_data_r[-drop,]
					VCV_r<-VCV_r[-drop, -drop]
				}
				pw_data_r<-droplevels(pw_data_r)
				
				# If we have the number of params missing work out
				if(is.na(other_params[[r]]) == T){
					other_params[[r]]<-length(unique(pw_data_r[,other_formulas[[r]]]))
				}
				
				# Test that there are sufficient data and that we do not just 1 level
				test_n<-((other_params[[r]] * 10) <= dim(pw_data_r)[1]) & (other_params[[r]] > 1)
				
				# run the meta-regression if there are enough data - otherwise just record NAs
				if(test_n == T){
					
					# for a numeric predictor fit a regular regression
					form<-paste0("~", other_formulas[[r]])
					# if the predictor is categorical remove the intercept
					if(is.numeric(pw_data_r[,other_formulas[[r]]]) == F){
						form<-paste0(form, " - 1")	
					}
					
					# Run the model	
					model<-try(rma.mv(yi=yi, V=VCV_r, random=list(~1|Experimental_Unit, ~1|ES_ID), data=pw_data_r, mods=as.formula(form)), silent=T)
				
				# Otherwise set the model as a dataframe
				}else{
					model<-data.frame(1)
				}
				# Also if the model fails set as a dataframe
				if(class(model)[1] == "try-error"){
					model<-data.frame(1)
				}

				# Save the reuslts for the rth regression in the mth imputation	
				trait_t_MRs_2[[r]][[m]]<-model
				
				# Clean up a bit
				rm(model)
				
				# Update the progress
				progress<-progress + 1
				setTkProgressBar(pb, progress)
			}
			
			# Clean up a bit
			rm(RMA)
			rm(pb_mod)
			rm(VCV)
			rm(low.covs, high.covs)
			rm(pw_data)
			rm(data_m)	
		}
		
		# Get the average for the imputed values and a copy pairwise dataset for plotting (we will use the average of imputed SDs)
		if(length(missing) > 0){
			av_imp<-sqrt(apply(my_imp^2, 1, mean))
			data$SD[as.numeric(names(av_imp))]<-av_imp
		}
		pw_data<-my_pairwise(data=data)
		pw_data<-escalc(n1i=exp_n, n2i=n, m1i=exp_Mean, m2i=Mean, sd1i=exp_SD, sd2i=SD, measure="ROM", data=pw_data)
		pw_data$diff_P<-pw_data$exp_Diet_P_kJ_g - pw_data$Diet_P_kJ_g
		pw_data$diff_BCAA<-pw_data$exp_Diet_BCAA_kJ_g - pw_data$Diet_BCAA_kJ_g
		pw_data$diff_nonBCAA<-pw_data$exp_Diet_nonBCAA_kJ_g - pw_data$Diet_nonBCAA_kJ_g
		pw_data$diff_Iso<-pw_data$exp_Diet_Iso_kJ_g - pw_data$Diet_Iso_kJ_g
		pw_data$diff_Leu<-pw_data$exp_Diet_Leu_kJ_g - pw_data$Diet_Leu_kJ_g
		pw_data$diff_Val<-pw_data$exp_Diet_Val_kJ_g - pw_data$Diet_Val_kJ_g
		pw_data$diff_BCAA.Non_BCAA_kJ_g<-pw_data$exp_Diet_BCAA.Non_BCAA_kJ_g - pw_data$Diet_BCAA.Non_BCAA_kJ_g
		
		# Save the results for the trait for MA and MR
		if(t==1){
			
			# Make a template list
			template<-lapply(traits, as.list)
			
			# Data
			data_list<-template
			data_list[[t]]<-pw_data
			
			# TSVs
			TSV_list<-template
			TSV_list[[t]]<-mean(trait_t_TSVs$TSV)
			
			# Meta Analysis
			MA_pool_list<-template
			MA_pool_list[[t]]<-my_agg(trait_t_MAs)
			
			# ER
			ER_pool_list<-template
			ER_pool_list[[t]]<-my_agg(trait_t_ERs)
			
			# ER
			TF_pool_list<-template
			TF_pool_list[[t]]<-my_agg(trait_t_TFs)
									
			# Meta-regression (nutrients)
			MR_pool_list<-template
			MR_pool_list[[t]]<-lapply(trait_t_MRs, my_agg)
			
			# Meta-regression (other)
			MR2_pool_list<-template
			MR2_pool_list[[t]]<-lapply(trait_t_MRs_2, my_agg)
						
		}else{
			data_list[[t]]<-pw_data
			TSV_list[[t]]<-mean(trait_t_TSVs$TSV)
			MA_pool_list[[t]]<-my_agg(trait_t_MAs)
			ER_pool_list[[t]]<-my_agg(trait_t_ERs)
			TF_pool_list[[t]]<-my_agg(trait_t_TFs)
			MR_pool_list[[t]]<-lapply(trait_t_MRs, my_agg)
			MR2_pool_list[[t]]<-lapply(trait_t_MRs_2, my_agg)	
		}
		
		# Clean up a bit
		rm(pw_data)
		rm(data)	
		rm(trait_t_TSVs)
		rm(trait_t_MAs)
		rm(trait_t_ERs)
		rm(trait_t_TFs)
		rm(trait_t_MRs)
		rm(trait_t_MRs_2)
		if(length(missing) > 0){
			rm(av_imp)
			rm(my_imp)
		}
	}	

	# Parcel up all the info and return it
	out<-list()
	out[[1]]<-data_list
	out[[2]]<-TSV_list
	out[[3]]<-MA_pool_list
	out[[4]]<-ER_pool_list
	out[[5]]<-TF_pool_list
	out[[6]]<-MR_pool_list
	out[[7]]<-MR2_pool_list
	names(out)<-c("Data_List", "TSV_List", "MA_List", "ER_List", "TF_List", "Nutri_MR_List", "Other_MR_List")
	close(pb)
	return(out)
	
}


