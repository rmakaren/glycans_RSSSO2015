
normalizations = function(dat) {
	# in this setting, for total area normalization we need exp(data) because they were simulated form log data
	glycans = grep("^GP\\d+$|^x\\d+$", names(dat), ignore.case = TRUE)
	
	expD = dat
	expD[, glycans] = exp(expD[, glycans])
	percD = totalAreaNorm(expD)
	logPercD = percD
	logPercD[,glycans] = log(logPercD[,glycans])
	logPercD$normalization = "logTA"
	
	RP = refPeakNorm(dat)
	Med = medianNorm(dat)
	MQ = medianQuotientNorm(dat)
	QN = quantileNorm(dat)
	
	percQN = quantileNorm(percD)
	percQN$normalization = "percQN"
	percMQ = medianQuotientNorm(percD)
	percMQ$normalization = "percMQ"
	rawQN = quantileNorm(expD)
	rawQN$normalization = "rawQN"
	
	dat$normalization = "Oracle"
	
	normData = data.frame(rbind(dat, percD, logPercD, RP, Med, MQ, QN, percQN, percMQ, rawQN))#
	return(normData)
}

totalAreaNorm = function(d) {	
	glycans = grep("^GP\\d+$|^x\\d+$", names(d), ignore.case = TRUE)
	tempD = d[, glycans]
	tempD = t(apply(tempD, 1,
                    function(samp){
                        samp/sum(samp)*100
                    }
                    ))

	d[,glycans] = tempD
	d$normalization = "TA"
	return(d)
}

refPeakNorm = function(f_data) {
	glycans = grep("^GP\\d+$|^x\\d+$", names(f_data), ignore.case = TRUE)
	d = f_data[glycans]
	max_peak = which.max(colSums(d, na.rm = TRUE)) 
	f_data[glycans] = round(d/d[,max_peak], digits = 5)
	f_data$normalization = "RP"
	return(f_data)
}

medianNorm = function(d) {
	glycans = grep("^GP\\d+$|^x\\d+$", names(d), ignore.case = TRUE)
	tempD = d[, glycans]
	tempD = t(apply(d[,glycans], 1,
                    function(gly){
                        (gly-median(gly, na.rm = TRUE))/IQR(gly, na.rm = TRUE)
                    }))
	
	d[,glycans] = tempD
	d$normalization = "Med"
	return(d)
}

medianQuotientNorm = function(f_data) {
	glycans = grep("^GP\\d+$|^x\\d+$", names(f_data), ignore.case = TRUE)
	ref_chrom = apply(f_data[glycans],2,median,na.rm = T)
	ref_chrom = sort(ref_chrom, decreasing = T)
	ref_chrom = as.matrix(f_data[,names(ref_chrom)]) %*% diag(1/ref_chrom)
	ref_chrom = apply(ref_chrom,1,median)
	f_data[,glycans] = round((as.matrix(f_data[,glycans]) / ref_chrom), digits = 2)
	
	f_data$normalization = "MQ"
	return(f_data)
}

quantileNorm = function(d) {
	glycans = grep("^GP\\d+$|^x\\d+$", names(d), ignore.case = TRUE)
	tempD = normalize.quantiles(t(as.matrix(d[, glycans])))
	
	d[,glycans] = t(tempD)
	d$normalization = "QN"
	return(d)
}


