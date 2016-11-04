#Script for spatial pattern analysis with wavelets and null models for transect data
#Two files required:
##One file with transect data - each column is a variable and each transect is a quadrat. The quadrats are contiguous
##One file with the sections (necessary for the CSRs, MC1s, AR2s simulations and for within-patch scale variance)


#Complete spatial randomness, no differences among sections - CSRh
CSRh <- function(x, Nperm=5000, keep.obs=T, print.loop=T, which.loop=100) {
	Nquad <- length(x)
	data.sim <- matrix(ncol=Nperm, nrow=Nquad)
	if(keep.obs) data.sim[,1] <- x
	for(i in ifelse(keep.obs, 2, 1):Nperm) {
		data.sim[,i] <- sample(x)
		if (print.loop) {
			if(i%%which.loop == 0) print(i)
		}
	}
	return(data.sim)
}

CSRs <- function(x, sections, Nperm=5000, keep.obs=T, print.loop=T, which.loop=100) {
	Nquad <- length(x)
	data.sim <- matrix(ncol=Nperm, nrow=Nquad)
	if(keep.obs) data.sim[,1] <- x
	for (i in ifelse(keep.obs, 2, 1):Nperm) {
		for (j in 1:nrow(sections)) {
			foo = sample(x[(sections$Start[j]:sections$End[j])])
			if (j == 1) bar = foo else bar = c(bar,foo)
			}
		data.sim[,i] <- bar
		if (print.loop) {
			if(i%%which.loop == 0) print(i)
		}
		
	}
	return(data.sim)
}

CSRd <- function(x, sections, disturbances=c("firebreak","firebreak_regenerating","railroad"), Nperm=5000, keep.obs=T, print.loop=T, which.loop=100) {
	Nquad <- length(x)
	data.sim <- matrix(ncol=Nperm, nrow=Nquad)
	if(keep.obs) data.sim[,1] <- x
	#Define which quadrats correspond to disturbed areas
	section.lengths <- sections$End-sections$Start+1
	data.section <- rep(as.character(sections[1,1]), section.lengths[1] )
	for (i in 2:nrow(sections)) {
		data.section <- c(data.section, rep(as.character(sections[i,1]), section.lengths[i]))
		}
	data.section <- as.factor(data.section)
	quad.disturb <- logical(Nquad)
	for (i in 1:Nquad) quad.disturb[i] <- any(data.section[i] == disturbances)
	quad.disturb <- which(quad.disturb)
	#Simulate the data
	for (i in ifelse(keep.obs, 2, 1):Nperm) {
		foo <- numeric(Nquad)
		foo[-quad.disturb] <- sample(x[-quad.disturb])
		foo[quad.disturb] <- sample(x[quad.disturb])
		data.sim[,i] <- foo
		if (print.loop) {
			if(i%%which.loop == 0) print(i)
		}
	}
	return(data.sim)
}

MC1h <- function(x, Nperm=5000, keep.obs=T, print.loop=T, which.loop=100) {
	Nquad <- length(x)
	data.sim <- matrix(ncol=Nperm, nrow=Nquad)
	if(keep.obs) data.sim[,1] <- x
	#Create transition matrix used in the MC1 simulations
	foo <- c(x[2:Nquad], NA)
	bar <- c(NA, x[1:(Nquad-1)])
	x.trans <- data.frame(x, foo, bar)
	foo <- c(x.trans$x, x.trans$x)
	bar <- c(x.trans$foo, x.trans$bar)	
	transmatrix <- table(foo,bar)
	sums <- apply(transmatrix, 1, sum)
	transmatrix <- transmatrix / sums
	for (i in ifelse(keep.obs, 2, 1):Nperm) {
		foo <- numeric(Nquad)
		foo.section <- x.trans$x
		foo.trans <- transmatrix
		foo.classes <- colnames(foo.trans)
		first <- sample(length(foo),1)
		foo[first] <- sample(foo.section,1)
		#Going backward
		if (first != 1 ) {			
			for(k in first:2) {
				ref <- as.character(foo[k])
				trans.row <- foo.trans[ref,]
				trans.row <- cumsum(trans.row)
				random <- runif(1)
				test <- random > trans.row
				foo[k-1] <- foo.classes[sum(test)+1]
			}
		}
		if (first != Nquad) {
			for(k in first:(length(foo)-1)) {
				ref <- as.character(foo[k])
				trans.row <- foo.trans[ref,]
				trans.row <- cumsum(trans.row)
				random <- runif(1)
				test <- random > trans.row
				foo[k+1] <- foo.classes[sum(test)+1]
			}
		}
		foobar <- as.numeric(foo)
		if (print.loop) {
			if(i%%which.loop == 0) print(i)
		}
		data.sim[,i] <- foobar
	}
	return(data.sim)
}

MC1s <- function(x, sections, Nperm=5000, keep.obs=T, print.loop=T, which.loop=100) {
	Nquad <- length(x)
	data.sim <- matrix(ncol=Nperm, nrow=Nquad)
	if(keep.obs) data.sim[,1] <- x
	section.lengths <- sections$End-sections$Start+1
	data.section <- rep(as.character(sections[1,1]), section.lengths[1])
	for (i in 2:nrow(sections)) {
		data.section <- c(data.section, rep(as.character(sections[i,1]), section.lengths[i]))
		}
	data.section <- as.factor(data.section)
	data.section.char <- as.character(data.section)
	foo <- c(x[2:Nquad], NA)
	bar <- c(NA, x[1:(Nquad-1)])
	x.trans <- data.frame(x, foo, bar)
	Nsections <- nrow(sections)	
	sectiontypes <- unique(data.section.char)
	Nsectiontypes <- length(sectiontypes)
	x.trans$loc <- character(Nquad)
	data.section.char.next <- c(data.section.char[2:Nquad],NA)
	data.section.char.prev <- c(NA, data.section.char[1:Nquad-1])
	#State whether is quadrat is located in the middle of a transect section or between sections; in this latter case it may not be included in the transition matrix.
	x.trans$loc <- ifelse(data.section.char == data.section.char.next & data.section.char == data.section.char.prev, "middle", ifelse(data.section.char != data.section.char.prev, "first", ifelse(data.section.char != data.section.char.next, "last", "error")))
	x.trans$loc[1] <- "first"
	x.trans$loc[Nquad] <- "last"
	#Make a transition matrix per section type
	trans.matrices <- list()
	for(i in 1:Nsectiontypes) {
		sectiontype <- sectiontypes[i]
		foobar.rows <- (1:Nquad)[data.section.char == sectiontype]
		foobar <- x.trans[foobar.rows,]
		foobar [foobar$loc == "last", 2] <- NA
		foobar[foobar$loc == "first", 3] <- NA
		foo <- c(foobar$x, foobar$x)
		bar <- c(foobar$foo, foobar$bar)	
		trans.matrices[[i]] <- table(foo,bar)
		sums <- apply(trans.matrices[[i]], 1, sum)
		trans.matrices[[i]] <- trans.matrices[[i]] / sums
	}		
	names(trans.matrices) = sectiontypes
	#Simulate the data, separately for each section
	for (i in ifelse(keep.obs, 2, 1):Nperm) {
		foobar = numeric()
		for (j in 1:Nsections) {
			foo.section <- as.character(sections[j,1])
			foo.trans <- trans.matrices[[foo.section]]
			foo.classes <- colnames(foo.trans)
			foo <- numeric((sections[j,3]-sections[j,2]+1) )
			first <- sample(length(foo),1)
			foo[first] <- sample(x[sections[j,2]:sections[j,3]],1)
			#Going backward
				if (first != 1 ) {			
					for(k in first:2) {
						ref <- as.character(foo[k])
						trans.row <- foo.trans[ref,]
						trans.row <- cumsum(trans.row)
						random <- runif(1)
						test <- random > trans.row
						foo[k-1] <- foo.classes[sum(test)+1]
					}
				}
				if (first != length(foo)) {
					for(k in first:(length(foo)-1)) {
						ref <- as.character(foo[k])
						trans.row <- foo.trans[ref,]
						trans.row <- cumsum(trans.row)
						random <- runif(1)
						test <- random > trans.row
						foo[k+1] <- foo.classes[sum(test)+1]
					}
				}
				if (j == 1) foobar = foo else foobar = c(foobar,foo)
			}
		data.sim[,i]=as.numeric(foobar)
		if (print.loop) {
			if(i%%which.loop == 0) print(i)
		}
	}
	return(data.sim)
}


MC1d <- function(x, sections, disturbances=c("firebreak","firebreak_regenerating","railroad"), Nperm=5000, keep.obs=T, print.loop=T, which.loop=100) {
	Nquad <- length(x)
	data.sim <- matrix(ncol=Nperm, nrow=Nquad)
	if(keep.obs) data.sim[,1] <- x
	#Define which quadrats correspond to disturbed areas
	section.lengths <- sections$End-sections$Start+1
	data.section <- rep(as.character(sections[1,1]), section.lengths[1] )
	for (i in 2:nrow(sections)) {
		data.section <- c(data.section, rep(as.character(sections[i,1]), section.lengths[i]))
		}
	quad.disturb <- logical(Nquad)
	for (i in 1:Nquad) quad.disturb[i] <- any(data.section[i] == disturbances)
	quad.disturb <- which(quad.disturb)
	#Classify quadrats as either disturbed or non-disturbed
	data.section[quad.disturb] <- "disturbance"
	data.section[-quad.disturb] <- "vegetation"
	data.section.next <- c(data.section[2:Nquad],NA)
	data.section.prev <- c(NA, data.section[1:Nquad-1])
	#Classify quadrats according to which part of the transect they're at; quadrats between disturbed areas and vegetation may not be included in the MC transition matrix.
	foo <- c(x[2:Nquad], NA)
	bar <- c(NA, x[1:(Nquad-1)])
	x.trans <- data.frame(x, foo, bar)
	x.trans$loc <- character(Nquad)
	x.trans$loc <- ifelse(data.section == data.section.next & data.section == data.section.prev, "middle", ifelse(data.section != data.section.prev, "first", ifelse(data.section != data.section.next, "last", "error")))
	x.trans$loc[1] <- "first"
	x.trans$loc[Nquad] <- "last"
	#Make transition matrix - consider all vegetation and homogeneous except for disturbances
	foobar <- x.trans[-quad.disturb,]
	foobar [foobar$loc == "last", 2] <- NA
	foobar[foobar$loc == "first", 3] <- NA
	foo <- c(x.trans$x, x.trans$x)
	bar <- c(x.trans$foo, x.trans$bar)	
	transmatrix <- table(foo,bar)
	sums <- apply(transmatrix, 1, sum)
	transmatrix <- transmatrix / sums
	#Simulate the data - first simulate an entire transect as in MC1h, then add disturbances as CSRd.
	for (i in ifelse(keep.obs, 2, 1):Nperm) {
		#The code here is longer than needed to maintain similarity with the MC1s code.
		foo <- numeric(Nquad)
		foo.section <- x.trans$x
		foo.trans <- transmatrix
		foo.classes <- colnames(foo.trans)
		first <- sample(length(foo),1)
		foo[first] <- sample(foo.section,1)
		#Going backward
		if (first != 1 ) {			
			for(k in first:2) {
				ref = as.character(foo[k])
				trans.row = foo.trans[ref,]
				trans.row = cumsum(trans.row)
				random = runif(1)
				test = random > trans.row
				foo[k-1] = foo.classes[sum(test)+1]
				}
			}
		if (first != length(foo)) {
			for(k in first:(length(foo)-1)) {
				ref <- as.character(foo[k])
				trans.row <- foo.trans[ref,]
				trans.row <- cumsum(trans.row)
				random <- runif(1)
				test <- random > trans.row
				foo[k+1] <- foo.classes[sum(test)+1]
				}
			}
		foobar <- as.numeric(foo)
		foobar[quad.disturb] <- sample(x[quad.disturb])
		data.sim[,i]=as.numeric(foobar)
		if (print.loop) {
			if(i%%which.loop == 0) print(i)
		}
	}
	return(data.sim)
}





