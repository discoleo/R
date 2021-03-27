########################
###
### Leonard Mada
### [the one and only]
###
### Percolation
###
### draft v.0.3c

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
revcol = function(m) {
	m = m[,rev(seq(ncol(m)))];
	invisible(m);
}
clean.percol = function(m, val=0) {
	# removes the non-percolating paths;
	pids = unique(m[,ncol(m)]);
	pids = pids[pids > 0];
	ids  = unique(m);
	ids  = ids[(ids > 0) & ! (ids %in% pids)]
	m[m %in% ids] = val;
	invisible(m);
}
shuffle.colors = function(m) {
	m2 = array(0, dim(m))
	m2[m < 0] = m[m < 0]
	vals = unique(as.vector(m))
	vals = vals[vals > 0];
	vnew = sample(vals, length(vals));
	id = 1
	for(val in vals) {
		m2[m == val] = vnew[id];
		id = id + 1;
	}
	invisible(m2)
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
	return(id)
}

norm.flux = function(m, add=0) {
	m[m > 0] = abs(log(add + m[m > 0]))
	invisible(m)
}

### Generators

rugrid.gen = function(dims, p) {
	m = sample(c(-1, 0), prod(dims), replace=T, prob=c(p, 1-p))
	m = matrix(m, nrow=dims[1])
}

rblock.gen = function(n, block.dim, min=0, max, prob, val=-1) {
	if(missing(max)) {
		if(length(min) >= 2) {
			max = min[2]; min = min[1];
		} else {
			max = prod(block.dim);
		}
	}
	val.count = seq(min, max);
	if(missing(prob)) {
		nn = sample(val.count, prod(n), replace=TRUE);
	} else {
		nn = sample(val.count, prod(n), replace=TRUE, prob=prob);
	}
	blocks.n = prod(block.dim);
	sample.block = function(n) {
		bm = array(as.integer(0), block.dim)
		npos = sample(seq(blocks.n), n)
		bm[npos] = val;
		bm
	}
	m = lapply(nn, sample.block);
	m = do.call(rbind, m);
	dims = n * block.dim;
	dim(m) = dims;
	invisible(m);
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
	invisible(m);
}

flood.all = function(m, type="Col1", val0=0, id.start, debug=TRUE) {
	type = match(type, c("Col1", "All"))
	if(is.na(type)) stop("Type NOT supported!")
	ncols = if(type == 2) seq(ncol(m)) else 1;
	
	if(missing(id.start)) {
		it = 1 + if(type == 2) max(m) else max(m[,1]);
	} else it = id.start;
	for(nc in ncols) {
		if(debug) print("Iteration: ");
		while(TRUE) {
			id.row = match(val0, m[,nc]);
			if(is.na(id.row)) break;
			if(debug) { cat(it); cat(", "); it = it + 1;}
			m = flood(m, c(id.row, nc), max(m)+1)
		}
		if(debug) cat("\n");
	}
	class(m) = c(class(m), if(type == 1) "filled" else "filledAll");
	invisible(m)
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
### Height
height.m = function(m) {
	ids = unique(as.vector(m));
	ids = ids[ids > 0];
	hm = array(as.integer(0), c(length(ids), 2, ncol(m)));
	
	for(nc in seq(ncol(m))) {
		for(nr in seq(nrow(m))) {
			val = m[nr, nc];
			val.id = match(val, ids);
			if(val <= 0) next;
			if(hm[val.id, 1, nc] == 0) {
				hm[val.id, 1, nc] = nr;
				hm[val.id, 2, nc] = nr;
			} else {
				hm[val.id, 2, nc] = nr;
			}
		}
	}
	invisible(hm);
}

### Count Area
count.fill = function(m, debug=TRUE) {
	isNotFilled = is.na(match("filledAll", class(m)))
	if(isNotFilled) m = flood.all(m, type="All");
	#
	count = table(m);
	ids   = as.integer(names(count));
	count = count[ids > 0];
	ids   = ids[ids > 0];
	if(debug) print(head(count), 20)
	# replace id with count:
	cn.m = m;
	for(id in seq(length(ids))) {
		cn.m[m == ids[id]] = count[id];
	}
	invisible(cn.m)
}

### Flow / Flux
flux = function(m, id, val0 = 1.0, debug=TRUE) {
	if(missing(id)) {
		id = max.id(m)
		if(debug) print(id);
	}
	#
	y.start = which(m[,1] %in% id)
	if(length(y.start) == 0) stop("NO such path!")
	vals = as.vector(rbind(y.start, 1, val0))
	# Init
	p.m = m;
	p.m[p.m != id] = -1;
	p.m[p.m == id] =  0;
	#
	pos = 1;
	# TODO: mixing of flows!
	while(pos <= length(vals)) {
		nn = double();
		while(pos <= length(vals)) {
			if(p.m[vals[pos], vals[pos + 1]] != 0) {pos = pos + 3; next;}
			p.m[vals[pos], vals[pos + 1]] = vals[pos + 2];
			fflow = 0;
			if(vals[pos] > 1 && p.m[vals[pos]-1, vals[pos + 1]] == 0) {
				nn = c(nn, vals[pos]-1, vals[pos + 1], 0);
				fflow = fflow + 1; }
			if(vals[pos] < nrow(m) && p.m[vals[pos]+1, vals[pos + 1]] == 0) {
				nn = c(nn, vals[pos]+1, vals[pos + 1], 0);
				fflow = fflow + 1; }
			if(vals[pos+1] > 1 && p.m[vals[pos], vals[pos + 1] - 1] == 0) {
				nn = c(nn, vals[pos], vals[pos + 1] - 1, 0);
				fflow = fflow + 1; }
			if(vals[pos+1] < ncol(m) && p.m[vals[pos], vals[pos + 1] + 1] == 0) {
				nn = c(nn, vals[pos], vals[pos + 1] + 1, 0);
				fflow = fflow + 1;
			}
			if(fflow == 0) {pos = pos + 3; next; }
			n.len = length(nn);
			nn[n.len - seq(0, fflow-1)*3] = vals[pos + 2] / fflow;
			pos = pos + 3;
		}
		vals = nn;
		pos = 1;
	}
	
	if(id != 0) p.m[m == 0] =  0;
	p.m[p.m < 0 & m > 0] =  0; # other non-connected "paths";
	return(p.m);
}
flux.dynamic = function(m, id, it=5, val0 = 1.0, debug=TRUE) {
	if(missing(id)) {
		id = max.id(m)
		if(debug) print(id);
	}
	#
	y.start = which(m[,1] %in% id)
	if(length(y.start) == 0) stop("NO such path!")
	vals = as.vector(rbind(y.start, 1, val0))
	# Init
	p.m = m;
	p.m[p.m != id] = -1;
	p.m[p.m == id] =  0;
	#
	pos = 1;
	# TODO: mixing of flows!
	for(itN in seq(it)) {
	while(pos <= length(vals)) {
		nn = double();
		while(pos <= length(vals)) {
			if(p.m[vals[pos], vals[pos + 1]] < 0) {pos = pos + 3; next;}
			# update Value
			valNew = p.m[vals[pos], vals[pos + 1]];
			valNew = valNew + vals[pos + 2];
			p.m[vals[pos], vals[pos + 1]] = valNew;
			cflow = 0; fflow = 0; # only push-forward Flow;
			if(vals[pos] > 1 && p.m[vals[pos]-1, vals[pos + 1]] >= 0) {
				valC = p.m[vals[pos]-1, vals[pos + 1]];
				if(valNew > valC + 1E-8) {
					nn = c(nn, vals[pos]-1, vals[pos + 1], valC);
					cflow = cflow + 1;
					fflow = fflow + valC;
				} }
			if(vals[pos] < nrow(m) && p.m[vals[pos]+1, vals[pos + 1]] >= 0) {
				valC = p.m[vals[pos]+1, vals[pos + 1]];
				if(valNew > valC + 1E-8) {
					nn = c(nn, vals[pos]+1, vals[pos + 1], valC);
					cflow = cflow + 1;
					fflow = fflow + valC;
				} }
			if(vals[pos+1] > 1 && p.m[vals[pos], vals[pos + 1] - 1] >= 0) {
				valC = p.m[vals[pos], vals[pos + 1] - 1];
				if(valNew > valC + 1E-8) {
					nn = c(nn, vals[pos], vals[pos + 1] - 1, valC);
					cflow = cflow + 1;
					fflow = fflow + valC;
				} }
			if(vals[pos+1] < ncol(m) && p.m[vals[pos], vals[pos + 1] + 1] >= 0) {
				valC = p.m[vals[pos], vals[pos + 1] + 1];
				if(valNew > valC + 1E-8) {
					nn = c(nn, vals[pos], vals[pos + 1] + 1, valC);
					cflow = cflow + 1;
					fflow = fflow + valC;
				}
			}
			if(cflow == 0) {pos = pos + 3; next; }
			# Update
			valNew = (valNew + fflow) / (cflow + 1);
			p.m[vals[pos], vals[pos + 1]] = valNew;
			n.len = length(nn); idNext = n.len - seq(0, cflow-1)*3;
			nn[idNext] = valNew - nn[idNext];
			pos = pos + 3;
		}
		vals = nn;
		pos = 1;
	}
		if(debug) print("Iteration")
		# TODO: evaluate NO new flow vs new flow;
		vals = as.vector(rbind(y.start, 1, 0));
	}
	
	if(id != 0) p.m[m == 0] =  0;
	p.m[p.m < 0 & m > 0] =  0; # other non-connected "paths";
	return(p.m);
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

### Flow
# - NO Mixing effects;
mflux = flux(m)
sum(mflux[m[,dim(m)[2]] > 0, dim(m)[2]])
mflux = norm.flux(mflux);
plot.rs(mflux)

# dynamic Flux:
# - with Mixing effects;
# - takes very LONG!!!
# TODO: debug!!!
mflux = flux.dynamic(m)
sum(mflux[m[,dim(m)[2]] > 0, dim(m)[2]])
mflux = norm.flux(mflux);
plot.rs(mflux)


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

plot.rs(split.rs(shuffle.colors(clean.percol(m))), main="Percolating Paths")

plot.rs(split.rs(clean.percol(m)), main="Percolating Paths")


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


### Average Path Height
# - may be useful for Normalization;
h.m = height.m(clean.percol(m))
h.diff.m = h.m[,2,] - h.m[,1,]
plot.rs(h.diff.m, "Height per column")
h.diff.m[1:10, 1:10]
apply(h.diff.m, 1, mean)


### Other:
# - Average Inputs;
# - Average Outputs;
# - Average Length;
# - Flux, Bottleneck;


path.all = length.path(reset.m(m), id=0)
plot.rs(split.rs(path.all), main="All Path Lengths")

### Mean Distance from In
# - includes also non-percolating paths;
mean(path.all[path.all > 0]) / ncol(path.all)
# ~ 0.98; # sum over the entire path length;
path.percol = path.all; path.percol[clean.percol(m) == 0] = 0;
mean(path.percol[path.percol > 0]) / ncol(path.percol)
# ~ 1.096; # sum over the entire path length;
### Mean Out Distance
nc = ncol(path.all)
mean(path.all[path.all[,nc] > 0, nc]) / nc
# ~ 2.02;


#####################
#####################

### Random Blocks
# - NOT uniformly random;

m = rblock.gen(c(450, 10), c(3,3))

plot.rs(split.rs(m))


path.m = flood.all(m)
path.m = shuffle.colors(path.m)
plot.rs(split.rs(path.m), main="Paths")



#############

### Ex 2
m = rblock.gen(c(450, 10), c(3,3), min=c(1,7))

plot.rs(split.rs(m))


path.m = flood.all(m)
path.m = shuffle.colors(path.m)
plot.rs(split.rs(path.m), main="Paths")



#############

### Ex 3:
# - foamy structure, little (or no) percolation;
m = rblock.gen(c(450, 10), c(3,3), min=c(4,5), prob=c(0.65, 0.35))

plot.rs(split.rs(m))


path.m = flood.all(m)
path.m = shuffle.colors(path.m)
plot.rs(split.rs(path.m), main="Paths")

table(path.m[,dim(m)[2]])

cfill.m = count.fill(path.m)
plot.rs(split.rs(cfill.m), main="Confluent Areas")

### Filter
# - remove small areas;
cfill.m[(cfill.m > 0) & (cfill.m < 5)] = 0
plot.rs(split.rs(cfill.m), main="Confluent Areas")

# How many inputs & In-widths:
table(path.m[,1])
# Area of inputs:
table(cfill.m[,1])


#############

### Ex 4:
# - foamy structure, little percolation;
m = rblock.gen(c(450, 10), c(3,3), min=c(4,5), prob=c(0.704, 0.296))

plot.rs(split.rs(m))


path.m = flood.all(m)
path.m = shuffle.colors(path.m)
plot.rs(split.rs(path.m), main="Paths")

table(path.m[,dim(m)[2]])

cfill.m = count.fill(path.m)
plot.rs(split.rs(cfill.m), main="Confluent Areas")

### Filter
# - remove small areas;
cfill.m[(cfill.m > 0) & (cfill.m < 5)] = 0
plot.rs(split.rs(cfill.m), main="Confluent Areas")

# How many inputs & In-widths:
table(path.m[,1])
# Area of inputs:
table(cfill.m[,1])


##############
##############

### Test
m = flood(m, c(1,1), max(m)+1)
m[1:10, 1:10]

table(m == 0)


### Test Raster
rs.m = toRaster(path.m);
plot(rs.m)

