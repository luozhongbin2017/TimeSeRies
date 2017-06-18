library(wmtsa)

wavCWTvarCIs.bivar <- function(x, y, significance=T, sections=NULL, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=T, CI=T, CIquant=c(0.025, 0.975), effect.size=T, keep.all=F, print.loop=T, which.loop=1, save.cwt=F, file.title="my_beautiful_wavelet", wav.template="gaussian2", scale.per.section=T) {
	#This function does not save the original wavelet output, retaining only scale and position variance. Such functionality may be added in the future. For now, the keep.all,save.cwt and file.title arguments are not fully functional.
	#This may not work for wavelets other than "gaussian2" ("Mexican hat") and "haar" (Haar).
	if (scale.per.section) {
		section.lengths <- sections$End-sections$Start+1
		data.section <- rep(as.character(sections[1,1]), section.lengths[1])
		for (k in 2:nrow(sections)) {
			data.section <- c(data.section, rep(as.character(sections[k,1]), section.lengths[k]))
		}
		data.section <- as.factor(data.section)
		section.types <- unique(sections[,1])
		Nsection.types <- length(section.types)
	}
	data.main.x <- x
	data.main.y <- y
	Nperm <- ifelse(test,5,ncol(data.main.x))
	dimnames(data.main.x) <- NULL
	dimnames(data.main.y) <- NULL
	data.foo.x <- data.main.x[,1]
	Nquad <- length(data.foo.x)
	full.length <- Nquad
	half.length <- round(Nquad/2)
	# Wavelet for x
	if(zero.padding) {
		data.foo.x <- c(rep(0,half.length),data.foo.x,rep(0,half.length/2)) #adds zeroes to the extremities of the transect
	}
	wav.orig.x <- wavCWT(data.foo.x,scale.range=c(1,scale.max),wavelet=wav.template)
	scales <- attr(wav.orig,"scale")
	Nscales <- length(scales)
	#Zero-padding adds zeroes to the extremities of the transect, to a length equal to half of the transect
	if(zero.padding) {
		#Save the attributes of the original wavelet transform
		wav.scale <- attr(wav.orig.x,"scale")
		wav.time <- 1:full.length
		wav.wavelet <- attr(wav.orig.x,"wavelet")
		wav.series <- attr(wav.orig.x,"series")
		wav.sampling.interval <- attr(wav.orig.x,"sampling.interval")
		wav.series.name <- attr(wav.orig.x,"series.name")
		wav.n.sample <- length(wav.time)
		wav.n.scale <- attr(wav.orig.x,"n.scale")
		wav.filter.arg <- attr(wav.orig.x,"filter.arg")
		#Remove parts of the wavelet transform corresponding to the artifficially added zeroes
		wav.orig.x <- wav.orig.x[(half.length+1):(half.length+full.length),]
		#Return the attributes to the corrected wavelet transform
		class(wav.orig.x) <- "wavCWT"
		attr(wav.orig.x,"scale") <- wav.scale
		attr(wav.orig.x,"time") <- wav.time
		attr(wav.orig.x,"wavelet") <- wav.wavelet
		attr(wav.orig.x,"series") <- wav.series
		attr(wav.orig.x,"sampling.interval") <- wav.sampling.interval
		attr(wav.orig.x,"series.names") <- wav.series.name
		attr(wav.orig.x,"n.sample") <- wav.n.sample
		attr(wav.orig.x,"n.scale") <- wav.n.scale
		attr(wav.orig.x,"filter.arg") <- wav.filter.arg
	}
	#The Cone Of Influence - COI - is the part of the wavelet transform that is affected by the transect's extremities, i.e. by quadrats for which no data has been collected. Inference made with these numbers is unreliable. The area affected is greater at larger scales, and therefore the COI has the shape of a, well, cone. Here I used an approach similar to that used by the PASSaGE software, which removes two quadrats for each scale for the MH wavelet and 1 quadrat for the Haar wavelet, from each side of the transect.
	if(remove.COI) {
		for(k in 1:length(scales)) {
			index1 <- 1:(scales[k]*ifelse(wav.template=="gaussian2",2,1))
			index2 <- nrow(wav.orig.x)-index1+1
			indices <- c(index1,index2)
			wav.orig.x[indices,k] <- NA
		}
	}
	# Wavelet for y 
	if(zero.padding) {
		data.foo.y <- c(rep(0,half.length),data.foo.y,rep(0,half.length/2)) #adds zeroes to the extremities of the transect
	}
	wav.orig.y <- wavCWT(data.foo.y,scale.range=c(1,scale.max),wavelet=wav.template)
	scales <- attr(wav.orig.y,"scale") # The scales should be the same for both wavelets.
	Nscales <- length(scales)
	#Zero-padding adds zeroes to the extremities of the transect, to a length equal to half of the transect
	if(zero.padding) {
		#Save the attributes of the original wavelet transform
		wav.scale <- attr(wav.orig.y,"scale")
		wav.time <- 1:full.length
		wav.wavelet <- attr(wav.orig.y,"wavelet")
		wav.series <- attr(wav.orig.y,"series")
		wav.sampling.interval <- attr(wav.orig.y,"sampling.interval")
		wav.series.name <- attr(wav.orig.y,"series.name")
		wav.n.sample <- length(wav.time)
		wav.n.scale <- attr(wav.orig.y,"n.scale")
		wav.filter.arg <- attr(wav.orig.y,"filter.arg")
		#Remove parts of the wavelet transform corresponding to the artifficially added zeroes
		wav.orig.y <- wav.orig.y[(half.length+1):(half.length+full.length),]
		#Return the attributes to the corrected wavelet transform
		class(wav.orig.y) <- "wavCWT"
		attr(wav.orig.y,"scale") <- wav.scale
		attr(wav.orig.y,"time") <- wav.time
		attr(wav.orig.y,"wavelet") <- wav.wavelet
		attr(wav.orig.y,"series") <- wav.series
		attr(wav.orig.y,"sampling.interval") <- wav.sampling.interval
		attr(wav.orig.y,"series.names") <- wav.series.name
		attr(wav.orig.y,"n.sample") <- wav.n.sample
		attr(wav.orig.y,"n.scale") <- wav.n.scale
		attr(wav.orig.y,"filter.arg") <- wav.filter.arg
	}
	#The Cone Of Influence - COI - is the part of the wavelet transform that is affected by the transect's extremities, i.e. by quadrats for which no data has been collected. Inference made with these numbers is unreliable. The area affected is greater at larger scales, and therefore the COI has the shape of a, well, cone. Here I used an approach similar to that used by the PASSaGE software, which removes two quadrats for each scale for the MH wavelet and 1 quadrat for the Haar wavelet, from each side of the transect.
	if(remove.COI) {
		for(k in 1:length(scales)) {
			index1 <- 1:(scales[k]*ifelse(wav.template=="gaussian2",2,1))
			index2 <- nrow(wav.orig.y)-index1+1
			indices <- c(index1,index2)
			wav.orig.y[indices,k] <- NA
		}
	}

	wav.covar <- as.matrix(wav.orig.x) * as.matrix(wav.orig.y)
	if(scale.var) wav.bivar.scale <- apply(wav.covar,2,mean, na.rm=T)
	if(pos.var) wav.bivar.pos <- apply(wav.covar,1,mean, na.rm=T)
	if(scale.per.section) {
		scale.var.section <- list()
		for (k in 1:Nsection.types) {
			quadrats.use <- (1:Nquad)[data.section==section.types[k]]
			scale.var.section[[k]] <- apply(wav.covar[quadrats.use,],2,mean,na.rm=T)
		}
		names(scale.var.section) <- section.types
	}
	if(scale.var) {
		data.scale <- matrix(ncol=ncol(data.main.x),nrow=Nscales)
		data.scale[,1] <- wav.bivar.scale
	}
	if(pos.var){
		data.pos <- matrix(ncol=ncol(data.main.x),nrow=length(wav.bivar.pos))
		data.pos[,1] <- wav.bivar.pos
	}
	if(scale.per.section){
		data.scale.section <- list()
		for (k in 1:Nsection.types) {
			data.scale.section[[k]] <- matrix(ncol=ncol(data.main.x), nrow=Nscales)
			data.scale.section[[k]][,1] <- scale.var.section[[k]]
		}
		names(data.scale.section) <- section.types
	}
	#The data.scale.section object stores the scale variance for each scale (scales are rows).
	if(save.cwt) write.table(wav.orig.x, paste(file.title,"_x_orig.txt",sep=""), sep=" ", dec=".", row.names=F, col.names=F)
	if(save.cwt) write.table(wav.orig.y, paste(file.title,"_y_orig.txt",sep=""), sep=" ", dec=".", row.names=F, col.names=F)
	#Significance calculation based on the randomized data. Set significance=F if no randomized values for significance calculation are supplied.
	if(significance) {
		for(j in 2:Nperm) {
			#The analysis of the randomized data starts here
			foo.x <- data.main.x[,j]
			full.length <- length(foo.x)
			half.length <- round(length(foo.x)/2)
			if(zero.padding) foo.x <- c(rep(0,half.length),foo.x,rep(0,half.length/2))
			foo.x.wav <- wavCWT(foo.x, scale.range=c(1,scale.max), wavelet=wav.template)
			if(zero.padding) {
				wav.scale <- attr(foo.x.wav,"scale")
				wav.time <- 1:full.length
				wav.wavelet <- attr(foo.x.wav,"wavelet")
				wav.series <- attr(foo.x.wav,"series")
				wav.sampling.interval <- attr(foo.x.wav,"sampling.interval")
				wav.series.name <- attr(foo.x.wav,"series.name")
				wav.n.sample <- length(wav.time)
				wav.n.scale <- attr(foo.x.wav,"n.scale")
				wav.filter.arg <- attr(foo.x.wav,"filter.arg")
				#Remove the uninformative wavelet transform values
				foo.x.wav <- foo.x.wav[(half.length+1):(half.length+full.length),]
				#Return the attributes of the wavelet transform
				class(foo.x.wav) <- "wavCWT"
				attr(foo.x.wav,"scale") <- wav.scale
				attr(foo.x.wav,"time") <- wav.time
				attr(foo.x.wav,"wavelet") <- wav.wavelet
				attr(foo.x.wav,"series") <- wav.series
				attr(foo.x.wav,"sampling.interval") <- wav.sampling.interval
				attr(foo.x.wav,"series.names") <- wav.series.name
				attr(foo.x.wav,"n.sample") <- wav.n.sample
				attr(foo.x.wav,"n.scale") <- wav.n.scale
				attr(foo.x.wav,"filter.arg") <- wav.filter.arg
			}
			if(remove.COI) {
				for(k in 1:length(scales)) {
					index1 <- 1:(scales[k]*ifelse(wav.template=="gaussian2",2,1))
					index2 <- nrow(wav.orig)-index1+1
					indices <- c(index1,index2)
					foo.x.wav[indices,k] <- NA
				}
			}
			if(all(c(keep.all, save.cwt))) write.table(foo.wav,paste(file.title,"_x_",formatC(j,digits=3,format="d",flag=0),".txt",sep=""),sep=" ", dec=".", row.names=F, col.names=F)

			foo.y <- data.main.y[,j]
			full.length <- length(foo.y)
			half.length <- round(length(foo.y)/2)
			if(zero.padding) foo.y <- c(rep(0,half.length),foo.y,rep(0,half.length/2))
			foo.y.wav <- wavCWT(foo.y, scale.range=c(1,scale.max), wavelet=wav.template)
			if(zero.padding) {
				wav.scale <- attr(foo.y.wav,"scale")
				wav.time <- 1:full.length
				wav.wavelet <- attr(foo.y.wav,"wavelet")
				wav.series <- attr(foo.y.wav,"series")
				wav.sampling.interval <- attr(foo.y.wav,"sampling.interval")
				wav.series.name <- attr(foo.y.wav,"series.name")
				wav.n.sample <- length(wav.time)
				wav.n.scale <- attr(foo.y.wav,"n.scale")
				wav.filter.arg <- attr(foo.y.wav,"filter.arg")
				#Remove the uninformative wavelet transform values
				foo.y.wav <- foo.y.wav[(half.length+1):(half.length+full.length),]
				#Return the attributes of the wavelet transform
				class(foo.y.wav) <- "wavCWT"
				attr(foo.y.wav,"scale") <- wav.scale
				attr(foo.y.wav,"time") <- wav.time
				attr(foo.y.wav,"wavelet") <- wav.wavelet
				attr(foo.y.wav,"series") <- wav.series
				attr(foo.y.wav,"sampling.interval") <- wav.sampling.interval
				attr(foo.y.wav,"series.names") <- wav.series.name
				attr(foo.y.wav,"n.sample") <- wav.n.sample
				attr(foo.y.wav,"n.scale") <- wav.n.scale
				attr(foo.y.wav,"filter.arg") <- wav.filter.arg
			}
			if(remove.COI) {
				for(k in 1:length(scales)) {
					index1 <- 1:(scales[k]*ifelse(wav.template=="gaussian2",2,1))
					index2 <- nrow(wav.orig)-index1+1
					indices <- c(index1,index2)
					foo.y.wav[indices,k] <- NA
				}
			}
			if(all(c(keep.all, save.cwt))) write.table(foo.wav,paste(file.title,"_y_",formatC(j,digits=3,format="d",flag=0),".txt",sep=""),sep=" ", dec=".", row.names=F, col.names=F)

			foo.covar <- as.matrix(foo.wav.x) * as.matrix(foo.wav.y)

			if(scale.var) foo.scale <- apply(foo.covar,2,mean,na.rm=T)
			if(pos.var) foo.pos <- apply(foo.covar,1,mean,na.rm=T)
			if(scale.per.section) {
				foo.var.section <- list()
				for (k in 1:Nsection.types) {
					quadrats.use <- (1:Nquad)[data.section==section.types[k]]
					foo.var.section[[k]] <- apply(as.matrix(foo.covar[quadrats.use,]),2,mean,na.rm=T)
				}
			}
			if(scale.var) data.scale[,j] <- foo.scale
			if(pos.var) data.pos[,j] <- foo.pos
			if(scale.per.section){
				for (k in 1:Nsection.types) {
					data.scale.section[[k]][,j] <- foo.var.section[[k]]
				}
			}
			if (print.loop) {
				if(j%%which.loop == 0) print(j)
			}
		}
		#This is the end of the loop. Each iteration calculated the wavelet transform on a randomized dataset and calculates scale covariance, position covariance, and scale covariance per section (if specified by the user).
		#Now calculation of confidence intervals and effect sizes.
		if (CI) { #Confidence intervals
			if(pos.var) pos.var.CI <- t(apply(data.pos, 1, quantile, CIquant, na.rm=T))
			if(scale.var) scale.var.CI <- t(apply(data.scale, 1, quantile, CIquant, na.rm=T))
			if(scale.per.section) {
				scale.section.CI <- list()
				for (k in 1:Nsection.types) {
					scale.section.CI[[k]] <- t(apply(data.scale.section[[k]], 1, quantile, CIquant, na.rm=T))
				}
				names(scale.section.CI) <- section.types
			}
		}
		if(effect.size) {
			if(pos.var) {
				pos.var.avg <- apply(data.pos, 1, mean, na.rm=T)
				pos.var.sd <- apply(data.pos, 1, sd, na.rm=T)
				pos.var.effect <- (wav.orig.pos - pos.var.avg) / pos.var.sd
			}
			if(scale.var) {
				scale.var.avg <- apply(data.scale, 1, mean, na.rm=T)
				scale.var.sd <- apply(data.scale, 1, sd, na.rm=T)
				scale.var.effect <- (wav.orig.scale - scale.var.avg) / scale.var.sd
			}
			if(scale.per.section) {
				scale.section.avg <- list()
				scale.section.sd <- list()
				scale.section.effect <- list()
				for (k in 1:Nsection.types) {
					scale.section.avg[[k]] <- apply(data.scale.section[[k]], 1, mean, na.rm=T)
					scale.section.sd[[k]] <- apply(data.scale.section[[k]], 1, sd, na.rm=T)
					scale.section.effect[[k]] <- (scale.var.section[[k]] - scale.section.avg[[k]]) / scale.section.sd[[k]]
				}
				names(scale.section.avg) <- section.types
				names(scale.section.sd) <- section.types
				names(scale.section.effect) <- section.types
			}
		}
		#Now, export the final object
		result <- list()
		k <- 1
		if (pos.var) {
			result[[k]] <- wav.orig.pos
			names(result)[k] <- "PosVar_Obs"
			k <- k+1
			if(CI) {
				result[[k]] <- pos.var.CI
				names(result)[k] <- "PosVar_CI" 
				k <- k+1
			}
			if(effect.size) { 
				result[[k]] <- data.frame(pos.var.avg, pos.var.sd, pos.var.effect)
				names(result[[k]]) <- c("Mean", "SD", "Effect_size")
				names(result)[k] <- "PosVar_Effect"
				k <- k+1
			}
		}
		if(any(c(scale.var, scale.per.section))) { 
			result[[k]] <- scales
			names(result)[[k]] <- "Scales"
			k <- k+1
		}
		if (scale.var) {
			result[[k]] <- wav.orig.scale
			names(result)[k] <- "ScaleVar_Obs"
			k <- k+1
			if(CI) { 
				result[[k]] <- scale.var.CI
				names(result)[k] <- "ScaleVar_CI" 
				k <- k+1
			}
			if(effect.size) {
				result[[k]] <- list(scale.var.avg, scale.var.sd, scale.var.effect)
				names(result[[k]]) <- c("Mean", "SD", "Effect_size")
				names(result)[k] <- "ScaleVar_Effect"
				k <- k+1
			}
		}
		if (scale.per.section) {
			result[[k]] <- scale.var.section
			names(result)[k] <- "ScaleSection_Obs"
			k <- k+1
			if(CI) {
				result[[k]] <-scale.section.CI 
				names(result)[k] <- "ScaleSection_CI" 
				k <- k+1
			}
			if(effect.size) {
				result[[k]] <- list(scale.section.avg, scale.section.sd, scale.section.effect)
				names(result[[k]]) <- c("Mean", "SD", "Effect_size")
				names(result)[k] <- "ScaleSection_Effect"
				k <- k+1
			}
		}
	} else {
		result <- list()
		k <- 1
		if(pos.var) {
			result[[k]] <- wav.orig.pos
			names(result)[k] <- "PosVar_Obs"
			k <- k+1
		}
		if (scale.var) {
		result[[k]] <- wav.orig.scale
			names(result)[k] <- "ScaleVar_Obs" 
			k <- k+1
		}
		if (scale.per.section) {
			result[[k]] <- scale.var.section
			names(result)[k] <- "ScaleSection_Obs"
			k <- k+1
		}
	}
	return(result)
}
