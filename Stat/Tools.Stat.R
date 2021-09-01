########################
###
### Leonard Mada
### [the one and only]
###
### Data Tools
###
### draft v.0.1a


### Tools for Basic Statistics


### Note:
# - some of the tools were created while mentoring students
#   and were used in their theses;
# - these are modified & improved versions;


######################

### Median / Quantiles

median.wt = function(x, wt, q=0.5, isSorted=FALSE) {
	# very simple implementation;
	# TODO:
	# 1.) https://aakinshin.net/posts/weighted-quantiles/
	# 2.) any(wt == 0)
	if(length(q) == 1) {
		if(q == 0) return(min(x[wt != 0]));
		if(q == 1) return(max(x[wt != 0]));
		if(q < 0 || q > 1) stop("Invalid quantile!");
	}
	if(any(q < 0 | q > 1)) stop("Invalid quantile!");
	# Sort
	if( ! isSorted) {
		id = order(x);
		x = x[id]; wt = wt[id];
	}
	# Median / Quantile
	wtcs = cumsum(wt);
	if(tail(wtcs, 1) == 0) stop("All weights are 0!");
	# Binary search
	search.bin = function(wtcs, cutoff, nstart=1, nend = length(wtcs)) {
		while(TRUE) {
			id = (nstart + nend) %/% 2;
			if(nend - nstart < 2) break; # TODO: exact
			if(wtcs[id] < cutoff) {
				nstart = id; # count = count + 1;
			} else if(wtcs[id] > cutoff) {
				nend = id; # count = count + 1;
			} else {
				break;
			}
		}
		# print(paste0("Count = ", count));
		return(id)
	}
	total = tail(wtcs, 1);
	# count = 0; # benchmarking
	if(length(q) > 1) {
		q.all.sorted = sort(unique(q));
		q.sorted = q.all.sorted;
		# Min/Max
		hasMin = if(q.sorted[1] == 0) {q.sorted = q.sorted[-1]; TRUE;} else FALSE;
		len = length(q.sorted);
		hasMax = if(q.sorted[len] == 1) {q.sorted = q.sorted[ - len]; TRUE;} else FALSE;
		# Quantiles
		id = sapply(q.sorted, function(q) search.bin(wtcs, total * q));
		# handle wt = 0 also in Min/Max;
		if(hasMin) {
			first = if(wtcs[1] != 0) 1 else
				for(i in seq(2, length(wtcs))) if(wtcs[i] != 0) return(i);
			id = c(first, id);
		}
		if(hasMax) {
			len = length(wtcs);
			last = if(wtcs[len] != 0) len else
				for(i in seq(len - 1, 1)) if(wtcs[i] != 0) return(i);
			id = c(id, last);
		}
		id = id[match(q, q.all.sorted)];
	} else {
		id = search.bin(wtcs, total * q);
	}
	return(x[id]);
}

