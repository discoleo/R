########################
###
### Leonard Mada
### [the one and only]
###
### Percolation
###
### draft v.0.1j

### Percolation

# - some experiments in Percolation

### Github:
# https://github.com/discoleo/R/blob/master/Stat/Percolation.R


####################
####################

### helper Functions

reset.m = function(m, id, val=0) {
	if(missing(id)) {
		m[m > 0] = val;
	} else if(id > 0) {
		m[m == id] = val;
	} else {
		m[(m > 0) & (m != id)] = val;
	}
	invisible(m)
}
### Percolation Functions
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

flood.all = function(m, type="by Col 1", val0=0, debug=TRUE) {
	# TODO: type
	while(TRUE) {
		if(debug) print("Iteration");
		id.row = match(val0, m[,1])
		if(is.na(id.row)) break;
		m = flood(m, c(id.row,1), max(m)+1)
	}
	return(m)
}

max.id = function(m) {
	out.m = m[, ncol(m)]
	out = table(out.m[out.m > 0])
	if(dim(out) == 0) {
		print("NO percolation!");
		out = table(m[m > 0])
	}
	id = match(max(out), out);
	id = as.integer(names(out)[id])
}
### Path Length
length.path = function(m, id, debug=TRUE) {
	if(missing(id)) {
		id = max.id(m)
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
	
	if(id != 0) p.m[m == 0] =  0;
	p.m[p.m < 0 & m > 0] =  0; # other non-connected "paths";
	return(p.m);
}
### Contact Area
contact.area = function(m, id) {
	area = 0;
	for(nc in 1:ncol(m)) {
		for(nr in 1:nrow(m)) {
			if(m[nr, nc] != id) next;
			if(nr > 1 && m[nr-1,nc] != id) area = area + 1;
			if(nr < nrow(m) && m[nr+1,nc] != id) area = area + 1;
			if(nc > 1 && m[nr,nc-1] != id) area = area + 1;
			if(nc < ncol(m) && m[nr,nc+1] != id) area = area + 1;
		}
	}
	return(area);
}

### Raster
toRaster = function(m, showVal=0) {
	rs.m = array(0, c(dim(m), 3));
	if( ! is.na(showVal)) {
		isZero = (m == showVal);
		doShow = TRUE;
	} else {
		doShow = FALSE;
	}

	### R
	layer.m = m;
	layer.m[m < 0] = 0
	val.max = max(layer.m);
	if(val.max > 0) layer.m = layer.m / val.max;
	if(doShow) layer.m[isZero] = 1;
	rs.m[,,1] = layer.m;

	### G
	layer.m = 1 - layer.m;
	layer.m[m <= 0] = 0
	if(doShow) layer.m[isZero] = 1;
	rs.m[,,2] = layer.m

	### B
	if(doShow) {
		layer.m = array(0, dim(m));
		layer.m[isZero] = 1;
		rs.m[,,3] = layer.m
	}

	rs.m = as.raster(rs.m)
	return(rs.m);
}
plot.rs = function(m, main, mar, line=0.5) {
	if( ! missing(main) ) hasTitle = TRUE else hasTitle = FALSE;
	if(missing(mar)) mar = c(0,0, if(hasTitle) 2 else 0, 0) + 0.1;
	type = match(class(m), c("raster", "matrix"));
	if(all(is.na(type))) stop("Data NOT supported!")
	if(any(type == 2)) {
		m = toRaster(m);
	}
	old.par = par(mar=mar);
		plot(m);
		if(hasTitle) mtext(main, line=line)
	par(old.par);
	invisible();
}
split.rs = function(m, n=5, from=1, max.len=5, w=10) {
	nr.tot = round(nrow(m) / n);
	frg.tot = ceiling(nrow(m) / nr.tot);
	# TODO: nrow(m) %% nr.tot > 0;
	if(from == 0) from = 1;
	if(from > 0) {
		frg.to = min(frg.tot, from + max.len);
	} else {
		from = max(1, frg.tot + 1 + from - max.len);
		frg.to = min(frg.tot, from + max.len);
	}
	m0 = matrix(0, ncol=w, nrow=nr.tot);
	m2 =  array(0, c(nr.tot, 0))
	for(frg in from:frg.to) {
		r.start = (frg - 1) * nr.tot + 1;
		r.end   = r.start + nr.tot - 1;
		m2 = cbind(m2, cbind(m[r.start:r.end,], m0))
	}
	invisible(m2);
}


################
################

### Examples
dims = c(80, 80)
p = 0.3

m = sample(c(-1, 0), prod(dims), replace=T, prob=c(p, 1-p))
m = matrix(m, nrow=dims[1])
m[1:10, 1:10]
plot.rs(m, "Percolation")


### Flood Fill
m = flood.all(m)

m[1:10, 1:10]

table(m)
table(m[,dims[2]])

plot.rs(m, "Percolation")


### Shortest Path
path.m = length.path(m)

id = dim(path.m)[2];
path.m[1:10, seq(id - 10, id)]

table(path.m[,dims[2]])


### Raster
plot.rs(rs.m, main="Path Length")


### Stat/Percolation
# - cells not accessible;
sum(m == 0) / prod(dim(m))
# - total contact area
a0 = contact.area(m, -1)
# - liquid contact area
a1 = contact.area(m, max.id(m))
a0; a1; a1 / a0;



###################

### Ex 2:
dims = c(80, 80)
p = 0.4

m = sample(c(-1, 0), prod(dims), replace=T, prob=c(p, 1-p))
m = matrix(m, nrow=dims[1])
m[1:10, 1:10]

m = flood.all(m)

m[1:10, 1:10]

table(m)
table(m[,dims[2]])


plot.rs(m, main="Percolation: Multiple Paths")


### Shortest Path
path.m = length.path(m)

id = dim(path.m)[2];
path.m[1:10, seq(id - 10, id)]

table(path.m[,dims[2]])


### Raster
plot.rs(path.m, main="Path Length")
# Note:
# - only "dominant" path is visualized;

# plot.rs(length.path(m, id=5), main="Path Length")


### Stat/Percolation
# - cells not accessible;
sum(m == 0) / prod(dim(m))
# - total contact area
a0 = contact.area(m, -1)
# - liquid contact area
a1 = contact.area(m, max.id(m))
a0; a1; a1 / a0;


#############

### Ex 3:
dims = c(200, 200) # takes long!
p = 0.40

m = sample(c(-1, 0), prod(dims), replace=T, prob=c(p, 1-p))
m = matrix(m, nrow=dims[1])
m[1:10, 1:10]

m = flood.all(m)

m[1:10, 1:10]

table(m)
table(m[,dims[2]])


#################
#################

### Periodic Boundary Condition

# - takes ~20 s;
dims = c(2000, 80)
p = 0.4

m = sample(c(-1, 0), prod(dims), replace=T, prob=c(p, 1-p))
m = matrix(m, nrow=dims[1])
m[1:10, 1:10]

m = flood.all(m)

m[1:10, 1:10]

table(m)
table(m[,dims[2]])


plot.rs(split.rs(m), main="Percolation: Multiple Paths")

# plot.rs(m, main="Percolation: Multiple Paths")


### Shortest Path
path.m = length.path(m)

id = dim(path.m)[2];
path.m[1:10, seq(id - 10, id)]

table(path.m[,dims[2]])


### Raster
plot.rs(path.m, main="Path Length")
# Note:
# - only "dominant" path is visualized;

# plot.rs(length.path(m, id=5), main="Path Length")

####################
### Stat/Percolation

### Inaccessible
# - cells that are not accessible;
sum(m == 0) / prod(dim(m))
# ~ 17%;
### TODO: compute over all paths;
# - total contact area
a0 = contact.area(m, -1)
# - liquid contact area
a1 = contact.area(m, max.id(m))
a0; a1; a1 / a0;


### Other:
# - Average Path Height;
# - Average Inputs;
# - Average Outputs;
# - Average Length;
# - Flux, Bottleneck;


path.all = length.path(reset.m(m), id=0)
plot.rs(split.rs(path.all), main="All Path Lengths")

### Mean Distance from In
# - includes also non-precolating paths;
mean(path.all[path.all > 0]) / ncol(path.all)
# ~ 0.98;
### Mean Out Distance
nc = ncol(path.all)
mean(path.all[path.all[,nc] > 0, nc]) / nc
# ~ 2.02;



###########

### Test
m = flood(m, c(1,1), max(m)+1)
m[1:10, 1:10]

table(m == 0)


### Test Raster
rs.m = toRaster(path.m);
plot(rs.m)

