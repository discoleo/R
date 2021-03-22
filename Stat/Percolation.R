########################
###
### Leonard Mada
### [the one and only]
###
### Percolation
###
### draft v.0.1b

### Percolation

# - some experiments in Percolation


####################

### helper Functions

flood = function(m, pyx, val=1, val0=0) {
	vals = pyx; pos = 1;
	
	while(TRUE) {
		if(pos > length(vals)) return(m);
		if(m[vals[pos], vals[pos + 1]] != val0) {pos = pos + 2; next;}
		m[vals[pos], vals[pos + 1]] = val;
		if(vals[pos] > 1) vals = c(vals, vals[pos]-1, vals[pos + 1]);
		if(vals[pos] < nrow(m)) vals = c(vals, vals[pos]+1, vals[pos + 1]);
		if(vals[pos+1] > 1) vals = c(vals, vals[pos], vals[pos + 1] - 1);
		if(vals[pos+1] < ncol(m)) vals = c(vals, vals[pos], vals[pos + 1] + 1);
		pos = pos + 2;
	}
	return(m);
}

flood.all = function(m, type="by Col 1", val0=0) {
	# TODO: type
	while(TRUE) {
		print("Iteration")
		id.row = match(val0, m[,1])
		if(is.na(id.row)) break;
		m = flood(m, c(id.row,1), max(m)+1)
	}
	return(m)
}

length.path = function(m, id, debug=TRUE) {
	if(missing(id)) {
		out.m = m[, ncol(m)]
		out = table(out.m[out.m > 0])
		id = match(max(out), out);
		id = as.integer(names(out)[id])
		if(debug) print(id);
	}
	p.m = m;
	p.m[p.m != id] = -1;
	p.m[p.m == id] =  0;
	# TODO
	lvl = 1; pos = 1;
	# rep(c(y, x))
	vals = as.vector(rbind(seq(nrow(m)), 1));
	while(pos <= length(vals)) {
		nn = integer();
		while(pos <= length(vals)) {
			if(p.m[vals[pos], vals[pos + 1]] != 0) {pos = pos + 2; next;}
			p.m[vals[pos], vals[pos + 1]] = lvl;
			if(vals[pos] > 1) nn = c(nn, vals[pos]-1, vals[pos + 1]);
			if(vals[pos] < nrow(m)) nn = c(nn, vals[pos]+1, vals[pos + 1]);
			if(vals[pos+1] > 1) nn = c(nn, vals[pos], vals[pos + 1] - 1);
			if(vals[pos+1] < ncol(m)) nn = c(nn, vals[pos], vals[pos + 1] + 1);
			pos = pos + 2;
		}
		vals = nn;
		lvl = lvl + 1; pos = 1;
	}
	
	p.m[m == 0] =  0;
	return(p.m);
}


### Examples
dims = c(80, 80)
p = 0.3

m = sample(c(-1, 0), prod(dims), replace=T, prob=c(p, 1-p))
m = matrix(m, nrow=dims[1])
m[1:10, 1:10]

m = flood.all(m)

m[1:10, 1:10]

table(m)
table(m[,dims[2]])

### Shortest Path
path.m = length.path(m)

id = dim(path.m)[2];
path.m[1:10, seq(id - 10, id)]

table(path.m[,dims[2]])


### Ex 2:
dims = c(200, 200) # takes long!
p = 0.40

m = sample(c(-1, 0), prod(dims), replace=T, prob=c(p, 1-p))
m = matrix(m, nrow=dims[1])
m[1:10, 1:10]

m = flood.all(m)

m[1:10, 1:10]

table(m)
table(m[,dims[2]])


###########

### Test
m = flood(m, c(1,1), max(m)+1)
m[1:10, 1:10]

table(m == 0)


