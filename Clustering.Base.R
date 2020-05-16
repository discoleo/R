#######################
###
### Advanced Clustering
### k-Means Extension
###
### Leonard Mada
###
### draft v 0.2

# Original on:
# https://github.com/discoleo/R/blob/master/Clustering.Base.R


# v 0.2
# - added cluster.apply function;
# - some experimentation with optimization;
# v 0.1b
# - Factors: are handled as factors;
# - Column Weights;
# - Metric: power specified as parameter;


######################
###  k-Means       ###
###  with Factors  ###
######################

### Distance functions

# the inlined version seems slightly faster;
# [possible because the functions have than less parameters]

# Dist Numeric
# non-inlined version:
dist.num.gen = function(p) {
	# optimized versions (but no metric)
	if(p == 1) {
		d.num = function(x, x.centroids, idCol, col.wt=1) {
			tmp.d = col.wt * (abs(x[,idCol] - x.centroids))
			return(tmp.d)
		}
	} else if(p == 2) {
		d.num = function(x, x.centroids, idCol, col.wt=1) {
			tmp.d = x[,idCol] - x.centroids;
			tmp.d = col.wt * tmp.d * tmp.d;
			return(tmp.d)
		}
	} else if(p - round(p) == 0 && p %% 2 == 0) {
		d.num = function(x, x.centroids, idCol, col.wt=1) {
			tmp.d = col.wt * (x[,idCol] - x.centroids)^p
			return(tmp.d)
		}
	} else {
		d.num = function(x, x.centroids, idCol, col.wt=1) {
			tmp.d = col.wt * (abs(x[,idCol] - x.centroids))^p
			return(tmp.d)
		}
	}
	return(d.num)
}
# Dist Factors
dist.fact.gen = function(p) {
	# p is currently NOT used;
	d.fact = function(x, x.centroids, idCol, col.wt=1) {
		# simple equality: dominant Factor
		level.centroid = levels(x[,idCol])[x.centroids]
		# tmp.d = ifelse(x[,idCol] == x.centroids, 0, col.wt)
		# slightly faster
		tmp.d = rep(col.wt, nrow(x))
		tmp.d[x[,idCol] == level.centroid] = 0
		# NOT faster
		# tmp.d = rep(col.wt, length(x))
		# tmp.d[x == x.centroids] = 0
		return(tmp.d)
	}
	return(d.fact)
}

### Cluster function
cluster = function(x, k, x.centroids, col.wt, p=2) {
	if(missing(col.wt)) {
		col.wt = rep(1, ncol(x))
	}
	### Distance Matrix
	# 1 column for each cluster;
	# k + 1 = save also Cluster ID for each data-row;
	# TODO: evaluate Matrix vs DF;
	d.df = matrix(0, nrow=nrow(x), ncol=k + 1)
	
	### Distance functions
	# Dist Numeric
	d.num = dist.num.gen(p);
	# Dist Factors
	d.fact = dist.fact.gen(p);
	
	num.Cols  = (1:ncol(x))[sapply(x[1,], function(x) is.numeric(x))]
	fact.Cols = (1:ncol(x))[sapply(x[1,], function(x) is.factor(x))]
	
	x.time = c(0, 0, 0, 0)

	### compute Distances
	id.rows = 1:nrow(d.df)
	for(i in 1:iter) {
		# Numeric
		x.time[1] = x.time[1] + system.time(
		for(idCol in num.Cols) {
			wt = col.wt[idCol]
			for(id.k in 1:k) {
				d.df[ , id.k] = d.df[ , id.k] + d.num(x, x.centroids[id.k, idCol], idCol, col.wt=wt)
			}
		}
		)[1]
		# seems: NO speed benefit
		# lapply(num.Cols,
			# function(idCol) {
				# wt = col.wt[idCol]
				# for(id.k in 1:k) {
					# d.df[,id.k] = d.df[,id.k] + d.num(x, x.centroids[id.k, idCol], idCol, col.wt=wt)
				# }
			# }
		# )
		# Factor
		x.time[2] = x.time[2] + system.time(
		for(idCol in fact.Cols) {
			wt = col.wt[idCol]
			for(id.k in 1:k) {
				d.df[ , id.k] = d.df[ , id.k] + d.fact(x, x.centroids[id.k, idCol], idCol, col.wt=wt)
			}
		}
		)[1]

		# summary(d.df)

		### assign Cluster
		d.df[ , k+1] = 1E+6 # HUGE VALUE
		d.df[ , k+1] = apply(d.df, M=1, min)
		new.id = sapply(id.rows, function(id) match(d.df[id, k+1], d.df[id, -(k+1)]))
		d.df[ , k+1] = new.id
	
		### new Centroids
		x.time[3] = x.time[3] + system.time(
		for(idCol in num.Cols) {
			for(id.k in 1:k) {
				# TODO: Mean, Median, ...
				x.centroids[id.k, idCol] = mean(x[ d.df[ , k+1] == id.k, idCol])
			}
		}
		)[1]
		#
		x.time[4] = x.time[4] + system.time(
		for(idCol in fact.Cols) {
			tmp.fact = levels(x.centroids[1, idCol])
			for(id.k in 1:k) {
				# TODO: calculate Factor Weights
				tmp = table(x[ d.df[ , k+1] == id.k, idCol])
				id.max = match(max(tmp), tmp)
				x.centroids[id.k, idCol] = tmp.fact[id.max]
			}
		}
		)[1]
		
		### reset Distances
		if(i == iter) {
			# last! # TODO
			cat("\n")
		} else {
			d.df[ - (k+1), ] = 0
			if(i %% 5 == 0) {
				cat(" "); cat(i);
				flush.console() # to display in real time
			}
		}
	}
	
	print(x.time) # Timings
	return(list("dist"=d.df, "centroids"=x.centroids))
}

### Cluster by Group function
cluster.apply = function(x, groups, k, x.centroids, col.wt, p=2) {
	if(missing(col.wt)) {
		col.wt = rep(1, ncol(x))
	}
	### Distance Matrix
	# 1 column for each cluster;
	# k + 1 = save also Cluster ID for each data-row;
	# TODO: evaluate Matrix vs DF;
	d.df = matrix(0, nrow=length(unique(groups)), ncol=k + 1)
	
	### Distance functions
	# Dist Numeric
	d.num = dist.num.gen(p);
	# Dist Factors
	d.fact = dist.fact.gen(p);
	
	num.Cols  = (1:ncol(x))[sapply(x[1,], function(x) is.numeric(x))]
	fact.Cols = (1:ncol(x))[sapply(x[1,], function(x) is.factor(x))]
	
	count.f = function(x) {
		sum( ! is.na(x) )
	}

	### compute Distances
	id.rows = 1:nrow(d.df)
	for(i in 1:iter) {
		# Numeric
		for(idCol in num.Cols) {
			wt = col.wt[idCol]
			for(id.k in 1:k) {
				tmp = d.num(x, x.centroids[id.k, idCol], idCol, col.wt=wt)
				d.df[,id.k] = d.df[,id.k] + tapply(tmp, groups, sum, na.rm=TRUE)
			}
		}
		# seems: NO speed benefit
		# lapply(num.Cols,
			# function(idCol) {
				# wt = col.wt[idCol]
				# for(id.k in 1:k) {
					# tmp = d.num(x, x.centroids[id.k, idCol], idCol, col.wt=wt)
					# d.df[,id.k] = d.df[,id.k] + tapply(tmp, groups, sum, na.rm=TRUE)
				# }
			# }
		# )
		# Factor
		for(idCol in fact.Cols) {
			wt = col.wt[idCol]
			for(id.k in 1:k) {
				tmp = d.fact(x, x.centroids[id.k, idCol], idCol, col.wt=wt)
				d.df[,id.k] = d.df[,id.k] + tapply(tmp, groups, sum, na.rm=TRUE)
			}
		}

		# summary(d.df)

		### assign Cluster
		d.df[, k+1] = 1E+6 # HUGE VALUE
		d.df[, k+1] = apply(d.df, M=1, min)
		new.id = sapply(id.rows, function(id) match(d.df[id, k+1], d.df[id, -(k+1)]))
		d.df[, k+1] = new.id
	
		### new Centroids
		# print("Updating Centroids")
		for(idCol in num.Cols) {
			# total.count = tapply(x[,idCol], d.df[groups, k+1], count.f)
			x.centroids[, idCol] = tapply(x[,idCol], d.df[groups, k+1], mean, na.rm=TRUE)
			# for(id.k in 1:k) {
				# TODO: Mean, Median, ...
				# x.centroids[id.k, idCol] = sum(x[ d.df[, k+1] == id.k, idCol])
			# }
		}
		# print("Finished Updating Numeric Centroids")
		for(idCol in fact.Cols) {
			tmp.fact = levels(x.centroids[1, idCol])
			for(id.k in 1:k) {
				# TODO: variant: compute Factor Weights
				tmp = table(x[ d.df[groups, k+1] == id.k, idCol])
				id.max = match(max(tmp), tmp)
				x.centroids[id.k, idCol] = tmp.fact[id.max]
			}
		}
		
		### reset Distances
		if(i == iter) {
			# last! # TODO
			cat("\n")
		} else {
			d.df[ , - (k+1)] = 0
			if(i %% 5 == 0) {
				cat(" "); cat(i);
				flush.console() # to display in real time
			}
		}
	}
	return(list("dist"=d.df, "centroids"=x.centroids))
}

