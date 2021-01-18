empirical_bayes_BC  <-  function(f_data, glycans=NULL) {
	if(is.null(glycans)){
		glycans = grep("^GP[0-9]+$",colnames(f_data),value=T)
	}
#	f_data_p = f_data["Age"]
	f_data_e = t(as.matrix(log(f_data[glycans])))
	f_data_b = f_data$Plate
	f_data_m = model.matrix(~ 1, data=f_data)
	combat_data = ComBat(dat=f_data_e, batch=f_data_b, mod=f_data_m,
			     par.prior=TRUE, prior.plots=FALSE) 
	f_data[glycans] = t(exp(combat_data))
	return(f_data)
}

