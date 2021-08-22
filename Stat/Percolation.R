########################
###
### Leonard Mada
### [the one and only]
###
### Percolation
###
### draft v.0.4a-cleanup

### Percolation

# - some experiments in Percolation
# - for Examples, see file:
#   Percolation.Examples.R;

### Github:
# https://github.com/discoleo/R/blob/master/Stat/Percolation.R


####################

###############
### History ###
###############


### draft v.0.4a-cleanup:
# - examples have been moved to file:
#   Percolation.Examples.R;


####################
####################

### helper Functions

# fast load
# source("Percolation.R")

tail.m = function(m, n=30, print=FALSE) {
	len = dim(m)[2];
	m = m[ , seq(max(1, len - n), len)]
	if(print) return(m);
	invisible(m);
}
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
rliniar.gen = function(n, w, d, ppore=3, pblock=0.5, val=-1) {
	# n = no. of channels; w = width;
	nc = d*n+n+1;
	m = matrix(0, nrow=w, ncol=nc);
	# Channel walls
	idChW = seq(1, nc, by=d+1)
	m[, idChW] = val;
	# add Pores
	npores = rpois(length(idChW), ppore)
	xwidth = seq(w);
	for(id in seq(n+1)) {
		xPore = sample(xwidth, npores[id])
		m[xPore, idChW[id]] = 0;
	}
	# block Channels
	nBlock = rpois(n, pblock)
	for(id in seq(n)) {
		if(nBlock[id] == 0) next;
		xBlock = sample(xwidth, nBlock[id]);
		# (id-1)*(d+1)+1
		m[xBlock, seq(idChW[id] + 1, length.out=d)] = -1;
	}
	
	return(t(m));
}
rlinwalk.gen = function(n, w, d, walk=c(-1,0,1), pwalk=c(1,2,1), ppore=3,
		first.const=TRUE, duplicate.at=NULL, val=-1) {
	# n = no. of channels; w = width of Material;
	# d = diameter of channel;
	nc = d*n+n+1;
	m = matrix(0, nrow=w, ncol=nc);
	# Channel walls:
	# walls start at fixed positions
	walk = walk * w + 1; # walk[walk == 0] = 1;
	nW = if(first.const) n else (n+1);
	wall = sample(walk, nW*(w-1), replace=TRUE, prob=pwalk);
	idChW = seq(1, nc, by=d+1);
	nposChWAbs = (idChW-1)*w + 1;
	wall = matrix(wall, nrow=w-1);
	if(first.const) wall = cbind(1, wall);
	wall = rbind(nposChWAbs, wall);
	wall = apply(wall, 2, cumsum);
	if( ! is.null(duplicate.at)) wall = cbind(wall, wall + duplicate.at*w);
	# unique(sort(...)): common paths may be rare;
	wall = unique(sort(as.vector(wall)));
	wall = wall[wall > 0 & wall <= (w*nc)];
	m[wall] = val;
	invisible(t(m));
	# TODO: add Pores
}
rrect.gen = function(n, dim, w.lim, h.lim, lambda.pores=2, addPores=TRUE,
		type=c("Poisson", "Simple"), aspect.fixed=NULL, prob.dir=c(1,1), val=-1) {
	# n = no. of rectangles;
	HD = 1; VD = 2;
	x0 = round(runif(n, 1, dim[HD]));
	y0 = round(runif(n, 1, dim[VD]));
	dw = round(runif(n, w.lim[1], w.lim[2]));
	if(is.null(aspect.fixed)) {
		dh = round(runif(n, h.lim[1], h.lim[2]));
	} else {
		dh = dw * aspect.fixed;
	}
	xdir.r = sample(c(-1,1), n, replace=TRUE, prob=prob.dir);
	ydir.r = sample(c(-1,1), n, replace=TRUE, prob=prob.dir);
	# OX
	xe = x0 + dw - 1; xs = x0 - dw; # offset: - 1;
	xe[xdir.r < 0] = x0[xdir.r < 0] - 1;
	xs[xdir.r > 0] = x0[xdir.r > 0] - 1;
	hasV2 = ifelse(xdir.r > 0, xe < dim[HD], xs >= 0);
	xe[xe >= dim[HD]] = dim[HD] - 1;
	xs[xs < 0] = 0;
	# OY
	ye = y0 + dh - 1; ys = y0 - dh + 1; # offset: -/+ 1 ???
	ye[ydir.r < 0] = y0[ydir.r < 0];
	ys[ydir.r > 0] = y0[ydir.r > 0];
	hasH2e = (ye <= dim[VD]);
	ye[ ! hasH2e] = dim[VD];
	hasH2s = (ys > 0);
	ys[ ! hasH2s] = 0;
	hasH2 = hasH2e & hasH2s;
	#
	vline.r = function(id) {
		l.seq = (ys[id]:ye[id]);
		if(hasV2[id]) {
			px = xs[id] * dim[VD] + l.seq;
			px = c(px, xe[id] * dim[VD] + l.seq);
		} else if(xdir.r[id] > 0) {
			px = xs[id] * dim[VD] + l.seq;
		} else {
			px = xe[id] * dim[VD] + l.seq;
		}
		return(px);
	}
	hline.r = function(id) {
		l.seq = (xs[id]:xe[id]) * dim[VD];
		if(hasH2[id]) {
			px = ys[id] + l.seq;
			px = c(px, ye[id] + l.seq);
		} else if(ydir.r[id] > 0) {
			px = ys[id] + l.seq;
		} else {
			px = ye[id] + l.seq;
		}
		return(px);
	}
	vl = unlist(lapply(seq(n), vline.r));
	hl = unlist(lapply(seq(n), hline.r));
	### Debug
	if(any(vl < 1)) print("ERROR: V");
	if(any(hl < 1)) print("ERROR: H");
	print(max(vl));
	print(max(hl));
	px = sort(unique(c(vl, hl)));
	m = matrix(0, nrow=dim[VD], ncol=dim[HD]);
	m[px] = val;
	### Pores
	if(addPores) {
		type = match.arg(type);
		if(type == "Simple") {
			nPores = round(n * lambda.pores);
			# nPores = sum(rpois(n, lambda.pores));
			print(paste0("Pores: ", nPores));
			idPores = sample(px, nPores);
			m[idPores] = 0;
		} else {
			pores = rpores(data.frame(hasH2 = hasH2, hasV2 = hasV2), lambda.pores);
		
			pxy = lapply(pores, function(x) {
				hasHV = attr(x, "f");
				xs = xs[x$id]; xe = xe[x$id]; ys = ys[x$id]; ye = ye[x$id];
				# x$Total: is NOT perfect;
				xm = round((xs * (x$Total + 1 - x$seq) + x$seq * xe)/ (x$Total + 1));
				ym = round((ys * (x$Total + 1 - x$seq) + x$seq * ye)/ (x$Total + 1));
				# xdir = xdir.r[x$id]; ydir = ydir.r[id];
				### Categories: OX = 1 & 2; OY = 3 & 4;
				if(all(hasHV[1,])) {
					print("Both")
					p = xm; # valid for both Cat: 1 & 2;
					p[x$cat == 3] = xs[x$cat == 3];
					p[x$cat == 4] = xe[x$cat == 4];
					p = p * dim[VD];
					p[x$cat >= 3] = p[x$cat >= 3] + ym[x$cat >= 3];
					p[x$cat == 1] = p[x$cat == 1] + ys[x$cat == 1];
					p[x$cat == 2] = p[x$cat == 2] + ye[x$cat == 2];
				} else if(hasHV$hasH2[1]) {
					print("H2")
					p = xm; # valid for both Cat: 1 & 2;
					p = p * dim[VD];
					p[x$cat == 1] = p[x$cat == 1] + ys[x$cat == 1];
					p[x$cat == 2] = p[x$cat == 2] + ye[x$cat == 2];
				} else if(hasHV$hasV2[1]) {
					p = ifelse(x$cat == 1, xs, xe);
					p = p * dim[VD] + ym;
				} else {
					# rectangles are in corners;
					# pores added to H1-stub;
					p = xm * dim[VD] + ifelse(ydir.r[x$id] > 0, ys, ye);
				}
				return(p);
			})
			pxy = sort(unique(unlist(pxy)));
			m[pxy] = 0; # set pores
			print(paste("Unique pores: ", length(pxy), collapse=""));
			return(m);
		### TODO: correct number of pores!
		}
	}
	invisible(m);
}
rpores = function(x, lambda) {
	n = nrow(x);
	nPores = rpois(n, lambda);
	hasPores = (nPores > 0);
	nPores = nPores[hasPores];
	### Categories
	# - using columns present in x;
	pores.df  = cbind(x[hasPores,], nPores = nPores);
	pores.cat = rpores.cat(pores.df);
	# by Category
	pores.df$id = seq(n)[hasPores];
	prle = split.cat(pores.df);
	# pl = list(pores=pores.df, hasPores=hasPores, cat=pores.cat, prle=prle);
	pl = lapply(seq(length(prle)), function(id) {
		x.df = prle[[id]];
		ids = rep(x.df$id, x.df$nPores);
		Total = rep(x.df$nPores, x.df$nPores);
		pores.df = data.frame(id=ids, cat = pores.cat$id.cat[[id]], Total=Total);
		p.seq = tapply.seq(pores.df[,c("id", "cat")]);
		pores.df$seq = p.seq;
		attr(pores.df, "f") = attr(x.df, "f");
		return(pores.df);
	})
	print(str(pl)) # Debug
	invisible(pl);
}
rpores.cat = function(pores.df) {
	# nPores per each Category:
	pores.total = aggregate(nPores ~ ., pores.df, sum);
	len = nrow(pores.total); # no. of categories;
	p.cat = data.frame(Any=TRUE, pores.total[ ! names(pores.total) %in% "nPores"]);
	p.cat.m = as.matrix(p.cat);
	id.cat = seq(ncol(p.cat.m)); # the categories;
	rp.cat = function(id) {
		sample(id.cat[p.cat.m[id,]], pores.total$nPores[id], replace=TRUE);
	}
	pores.l = list(p.cat=p.cat, id.cat=lapply(seq(len), rp.cat));
	return(pores.l);
}
split.cat = function(x, var.names=c("id", "nPores"), aggr.var=var.names[1]) {
	names.p = names(x);
	names.p = names.p[ ! names.p %in% var.names];
	# only in R 4.1.0!
	# frml = reformulate(names.p)
	prle = split(x[, var.names], f=x[ , names.p], drop=TRUE)
	var1 = var.names[1];
	aggr = aggregate(
		formula(paste0(aggr.var, "~.")),
		x[,c(aggr.var, names.p)], length)
	names(aggr)[length(names(aggr))] = "len";
	for(id in seq(length(prle))) {
		attr(prle[[id]], "f") = aggr[id,];
	}
	prle;
}
tapply.seq = function(x, id1=1, FUN=seq_along) {
	ave(x[,id1], x, FUN=FUN);
}

TEST = FALSE;
if(TEST) {
	### TODO: multiple pores by category;
	# m = rrect.gen(120, c(40, 200), c(6, 20), c(6, 16), lambda=3, aspect.fixed=2)
	m = rrect.gen(120, c(40, 200), c(6, 20), c(6, 16), lambda=3)
	plot.rs(m)

	m.fl = flood.all(m)
	plot.rs(m.fl)
}

#########################
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

### Diffusion
# - simple diffusion;
diffusion = function(m, id, val0 = 1.0, debug=TRUE) {
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
	p.m = diffusion.internal(p.m, vals)
	
	if(id != 0) p.m[m == 0] =  0;
	p.m[p.m < 0 & m > 0] =  0; # other non-connected "paths";
	invisible(p.m);
}
diffusion.internal = function(p.m, vals) {
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
	invisible(p.m);
}
### Dynamic diffusion [old]
diffusion.dynamic.slow = function(m, id, iter=5, val0 = 1.0, max.size.scale=3, debug=TRUE) {
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
	# pre-compute Diffusion
	# - better results, but still extremely slow!
	p.m = diffusion.internal(p.m, vals);
	# next iterations
	vals = as.vector(rbind(y.start, 1, 0)); # 0 vs val0?
	pos = 1; tol = 1E-24;
	for(itN in seq(iter)) {
	while(pos <= length(vals)) {
		nn = double();
		while(pos <= length(vals)) {
			if(p.m[vals[pos], vals[pos + 1]] < 0) {pos = pos + 3; next;}
			# update Value
			valNew = p.m[vals[pos], vals[pos + 1]];
			valNew = valNew + vals[pos + 2];
			p.m[vals[pos], vals[pos + 1]] = valNew;
			cflow = 0; fflow = 0; # only push-forward Flow;
			nnew = double();
			if(vals[pos+1] < ncol(m) && p.m[vals[pos], vals[pos + 1] + 1] >= 0) {
				valC = p.m[vals[pos], vals[pos + 1] + 1];
				if(valNew > valC + tol) {
					nnew = c(nnew, vals[pos], vals[pos + 1] + 1, valC);
					cflow = cflow + 1;
					fflow = fflow + valC;
				} }
			if(vals[pos] < nrow(m) && p.m[vals[pos]+1, vals[pos + 1]] >= 0) {
				valC = p.m[vals[pos]+1, vals[pos + 1]];
				if(valNew > valC + tol) {
					nnew = c(nnew, vals[pos]+1, vals[pos + 1], valC);
					cflow = cflow + 1;
					fflow = fflow + valC;
				} }
			if(vals[pos+1] > 1 && p.m[vals[pos], vals[pos + 1] - 1] >= 0) {
				valC = p.m[vals[pos], vals[pos + 1] - 1];
				if(valNew > valC + tol) {
					nnew = c(nnew, vals[pos], vals[pos + 1] - 1, valC);
					cflow = cflow + 1;
					fflow = fflow + valC;
				} }
			if(vals[pos] > 1 && p.m[vals[pos]-1, vals[pos + 1]] >= 0) {
				valC = p.m[vals[pos]-1, vals[pos + 1]];
				if(valNew > valC + tol) {
					nnew = c(nnew, vals[pos]-1, vals[pos + 1], valC);
					cflow = cflow + 1;
					fflow = fflow + valC;
				}
			}
			if(cflow == 0) {pos = pos + 3; next; }
			# Update
			valNew = (valNew + fflow) / (cflow + 1);
			n.len = length(nnew); idNext = n.len - seq(0, cflow-1)*3;
			# Priority: lower initial value;
			n.id = order(nnew[idNext]);
			idSorted = idNext[n.id];
			nnew = nnew[rbind(idSorted-2, idSorted-1, idSorted)];
			# Propagate Update immediately
			p.m[vals[pos], vals[pos + 1]] = valNew;
			for(idN1 in idNext) {
				p.m[nnew[idN1 - 2], nnew[idN1 - 1]] = valNew;
			}
			nnew[idNext] = 0;
			nn = c(nn, nnew);
			pos = pos + 3;
			if(length(nn) > max.size.scale * prod(dim(m))) {
				print("Internal Break!")
				break;
			}
		}
		vals = nn;
		pos = 1;
	}
		if(debug) print(paste0("Iteration: ", itN));
		# TODO: evaluate NO new flow vs new flow;
		vals = as.vector(rbind(y.start, 1, 0));
		pos = 1;
	}
	
	if(id != 0) p.m[m == 0] =  0;
	p.m[p.m < 0 & m > 0] =  0; # other non-connected "paths";
	return(p.m);
}
### Dynamic diffusion:
# - new Sequential Algorithm;
diffusion.dynamic = function(m, id, iter=40, val0 = 1.0, debug=TRUE) {
	if(missing(id)) {
		id = max.id(m)
		if(debug) print(id);
	}
	#
	y.start = which(m[,1] %in% id)
	if(length(y.start) == 0) stop("NO such path!")
	# Init
	p.m = m;
	p.m[m != id] = -1;
	p.m[m == id] =  0;
	p.m[y.start, 1] = val0; # start of flow;
	tol = 1E-24;
	#
	for(itN in seq(iter)) {
		for(nc in seq(ncol(m))) {
		for(nr in seq(nrow(m))) {
			if(p.m[nr, nc] <= 0) next;
			valNew = p.m[nr, nc];
			nn = double();
			cflow = 0; fflow = 0;
			if(nc < ncol(m) && p.m[nr, nc + 1] >= 0) {
				valC = p.m[nr, nc + 1];
				if(valNew > valC + tol) {
					nn = c(nn, nr, nc + 1);
					cflow = cflow + 1;
					fflow = fflow + valC;
				} }
			if(nr < nrow(m) && p.m[nr + 1, nc] >= 0) {
				valC = p.m[nr + 1, nc];
				if(valNew > valC + tol) {
					nn = c(nn, nr + 1, nc);
					cflow = cflow + 1;
					fflow = fflow + valC;
				} }
			if(nc > 1 && p.m[nr, nc - 1] >= 0) {
				valC = p.m[nr, nc - 1];
				if(valNew > valC + tol) {
					nn = c(nn, nr, nc - 1);
					cflow = cflow + 1;
					fflow = fflow + valC;
				} }
			if(nr > 1 && p.m[nr - 1, nc] >= 0) {
				valC = p.m[nr - 1, nc];
				if(valNew > valC + tol) {
					nn = c(nn, nr - 1, nc);
					cflow = cflow + 1;
					fflow = fflow + valC;
				}
			}
			if(cflow == 0) next;
			# Update
			valNew = (valNew + fflow) / (cflow + 1);
			p.m[nr, nc] = valNew;
			n.len = length(nn); idNext = n.len - seq(0, cflow-1)*2;
			for(idNext in (n.len - seq(0, cflow-1)*2)) {
				p.m[nn[idNext - 1], nn[idNext]] = valNew;
			}
		}
		}
		if(debug) print(paste0("Iteration: ", itN))
		# add new flow;
		p.m[y.start, 1] = p.m[y.start, 1] + val0;
	}
	
	if(id != 0) p.m[m == 0] =  0;
	p.m[p.m < 0 & m > 0] =  0; # other non-connected "paths";
	return(p.m);
}

#######################
### Graphical functions

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
	# w = width between displayed fragments;
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
	m0 = matrix(0, ncol=w, nrow=nr.tot); # spacer
	m2 = array(0, c(nr.tot, 0));
	for(frg in from:frg.to) {
		r.start = (frg - 1) * nr.tot + 1;
		r.end   = r.start + nr.tot - 1;
		if(r.end > nrow(m)) {
			m1 = matrix(0, ncol=ncol(m), nrow= r.end - nrow(m))
			m2 = cbind(m2, rbind(m[r.start:nrow(m),], m1), m0)
		} else {
			m2 = cbind(m2, m[r.start:r.end,], m0)
		}
	}
	invisible(m2);
}

###############
### Geo-Physics

### Bulk Modulus
mod.bulk = function(phi, K) {
	# Gassman's eq:
	# Ksat/(Kmin - Ksat) = Kdry/(Kmin - Kdry) + Kfl/(phi*(Kmin - Kfl))
	# (K[mineral], K[dry], K[fluid])
	if(length(K) < 3) {
		# Kfl: use Batzle & Wang's eqs;
	}
	dK = (K[1] - K[2:3]);
	Ks = K[2]/dK[2] + K[3]/(dK[3] * phi);
	Ksat = K[1] * Ks / (Ks + 1);
	return(Ksat);
}
mod.dry.bulk = function(phi, K) {
	# Kdry = (Ksat*(phi*Kmin/Kfl + 1 - phi) - Kmin) /
	#        (phi*Kmin/Kfl + Ksat/Kmin - 1 - phi);
	# Ksat[initial] = rho * (vp^2 - 4/3*vs^2);
}
mod.fluid.bulk = function(K, p) {
	# mixture of oil/fluids: p = proportion;
	if(length(K) - length(p) == 1) p = c(p, 1 - sum(p));
	1 / sum(p/K);
}
mod.min.bulk = function(type=c("quartz", "calcite", "dolomite", "muscovite",
		"feldspar", "albite", "halite", "anhydrite", "pyrite", "siderite")) {
	if(missing(type)) stop("Mineral name needed!")
	type = pmatch(type, formals(mod.min.bulk)$type[-1]);
	if(is.na(type)) stop("Mineral name not yet available!")
	
	K = c(36.6, 76.8, 94.9, 61.5, 75.6, 59.5, # albite: simulated 55;
		24.8, 56.1, 147.4, 123.7)
	return(K[type]);
	# albite:
	# https://www.ceramics-silikaty.cz/2015/pdf/2015_04_326.pdf
}


################
################

### Examples:
# - moved to file:
#   Percolation.Examples.R;
