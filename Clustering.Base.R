#######################
###
### Advanced Clustering
### k-Means Extension
###
### Leonard Mada
###
### draft v 0.1

# v 0.1
# - Factors: are handled as factors;
# - Weights;
# - Metric: power as parameter;


##########################
### k-Means on Factors ###
##########################

### Distance functions
# Dist Numeric
dist.num.gen = function(p) {
	if(p == 1) {
		d.num = function(x, d.df, x.centroids, idCol, centroid.id, wt=1) {
			tmp = x.centroids[centroid.id, idCol]
			tmp.d = d.df[ , centroid.id] + wt * (abs(x[,idCol] - tmp))
			return(tmp.d)
		}
	} else if(p == 2) {
		d.num = function(x, d.df, x.centroids, idCol, centroid.id, wt=1) {
			tmp = x.centroids[centroid.id, idCol]
			tmp.d = x[,idCol] - tmp;
			tmp.d = d.df[ , centroid.id] + wt * tmp.d * tmp.d;
			return(tmp.d)
		}
	} else if(p - round(p) == 0 && p %% 2 == 0) {
		d.num = function(x, d.df, x.centroids, idCol, centroid.id, wt=1) {
			tmp = x.centroids[centroid.id, idCol]
			tmp.d = d.df[ , centroid.id] + wt * (x[,idCol] - tmp)^p
			return(tmp.d)
		}
	} else {
		d.num = function(x, d.df, x.centroids, idCol, centroid.id, wt=1) {
			tmp = x.centroids[centroid.id, idCol]
			tmp.d = d.df[ , centroid.id] + wt * (abs(x[,idCol] - tmp))^p
			return(tmp.d)
		}
	}
	return(d.num)
}
# Dist Factors
dist.fact.gen = function(p) {
	# p is currently NOT used;
	d.fact = function(x, d.df, x.centroids, idCol, centroid.id, wt=1) {
		# simple equality: dominant Factor
		tmp = x.centroids[centroid.id, idCol]
		tmp.d = d.df[ , centroid.id] + ifelse(x[,idCol] == tmp, 0, wt)
		return(tmp.d)
	}
	return(d.fact)
}

### Cluster function
cluster = function(x, k, x.centroids, x.wt, p=2) {
	if(missing(x.wt)) {
		x.wt = rep(1, ncol(x))
	}
	### Distance Matrix
	# k + 1 = save also Cluster ID;
	d.df = matrix(0, nrow=nrow(x), ncol=k + 1)
	
	### Distance functions
	# Dist Numeric
	d.num = dist.num.gen(p);
	# Dist Factors
	d.fact = dist.fact.gen(p);
	
	num.Cols  = (1:ncol(x))[sapply(x[1,], function(x) is.numeric(x))]
	fact.Cols = (1:ncol(x))[sapply(x[1,], function(x) is.factor(x))]

	### compute Distances
	id.rows = 1:nrow(d.df)
	for(i in 1:iter) {
		# Numeric
		for(idCol in num.Cols) {
			wt = x.wt[idCol]
			for(id.k in 1:k) {
				d.df[,id.k] = d.num(x, d.df, x.centroids, idCol, id.k, wt=wt)
			}
		}
		# Factor
		for(idCol in fact.Cols) {
			wt = x.wt[idCol]
			for(id.k in 1:k) {
				d.df[,id.k] = d.fact(x, d.df, x.centroids, idCol, id.k, wt=wt)
			}
		}

		# summary(d.df)

		### assign Cluster
		d.df[, k+1] = 1E+6 # HUGE VALUE
		d.df[, k+1] = apply(d.df, M=1, min)
		new.id = sapply(id.rows, function(id) match(d.df[id, k+1], d.df[id, -(k+1)]))
		d.df[, k+1] = new.id
	
		### new Centroids
		for(idCol in num.Cols) {
			for(id.k in 1:k) {
				x.centroids[id.k, idCol] = mean(x[d.df[, k+1] == id.k, idCol])
			}
		}
		for(idCol in fact.Cols) {
			tmp.fact = levels(x.centroids[1, idCol])
			for(id.k in 1:k) {
				# TODO: with Factor Weights
				tmp = table(x[d.df[, k+1] == id.k, idCol])
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
	return(d.df)
}

