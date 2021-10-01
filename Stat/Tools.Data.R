########################
###
### Leonard Mada
### [the one and only]
###
### Data Tools
###
### draft v.0.1o-fix


### Tools to Process/Transform Data


###############
### History ###
###############


### draft v.0.1o - v.0.1o-fix:
# - basic implementation of a rename() function;
# - sequential renaming; [v.0.1o-fix]
### draft v.0.1n:
# - major refactoring of function ftable2;
### draft v.0.1l - v.0.1m:
# - [changed] justify: argument as list; [v.0.1l]
# - adapted expand.labels (makeLabels) function; [v.0.1m]
### draft v.0.1k [fix-7]:
# - [fix] align = center; [fix-1]
# - [fix] proper argument: split.ch; [fix-2]
# - [fix] added missing argument: pos; [fix-3]
# - [fix] added: method = row.compact & col.compact; [fix-4 & 5]
# - [fix] proper justifying; [fix-6]
# - [implemented] justify = center; [fix-7]
### draft v.0.1j:
# - the Section on Formulas/Expressions
#   has been moved to a separate file:
#   Tools.Formulas.R;


#######################
#######################

### Groups / Aggregates

# TODO


### DF
# - rename columns: simple implementation;
# - dplyr::rename: behaves differently during name clashes!
# Note:
# - implementations of the form (new.name = old.name):
#   do NOT permit the use of formulas to generate dynamically the new names;
rename2 = function(x, ...) {
	e = substitute(c(...));
	len = length(e);
	if(len <= 1) {
		print("No names were renamed!");
		return(x);
	}
	len = len - 1;
	e = e[-1];
	nms = names(e);
	isDuplicated = duplicated(nms);
	if(any(isDuplicated))
		stop(paste0("Duplicated names: ", paste0(nms[isDuplicated], collapse=", ")));
	# Old Names
	nms.old = character(len);
	for(id in seq(len)) {
		tmp.e = e[[id]];
		if( ! nzchar(nms[[id]])) {
			if(is.call(tmp.e)) nms[[id]] = as.character(eval(tmp.e[[2]]))
			else stop("Error: Missing name!")
		}
		if(is.call(tmp.e)) {
			tmp.e = eval(tmp.e);
			if(length(tmp.e) > 1) {
				stop("Multiple names mapped to same name!")
			}
		}
		print(tmp.e)
		nms.old[id] = as.character(tmp.e);
	}
	tmp = x; tmp.nms = names(tmp);
	# Sequential processing: as in sequential renaming;
	for(id in seq(length(nms))) {
		idn = match(nms.old[id], tmp.nms);
		if(is.na(idn)) stop(paste0("Error: Some of the names were NOT found!\n  Name: ",
			nms.old[id], " (possibly twice renamed)"));
		if(tmp.nms[idn] == nms[id]) next;
		if(any(tmp.nms == nms[id]))
			stop(paste0("Error: New names clash with the old names!\n  Names: ", nms[id]));
		tmp.nms[idn] = nms[id];
	}
	names(tmp) = tmp.nms;
	return(tmp);
}

### Row DF
# Matrix => Row DF
toRow.df = function(m, val=rownames(m), group=colnames(m), asNumVal=TRUE, asNumGroup=TRUE) {
	if(asNumVal && (! is.numeric(val))) val = as.numeric(val);
	if(asNumGroup && (! is.numeric(group))) group = as.numeric(group);
	count.df = data.frame(
		v = val,
		count = as.vector(m),
		group = rep(group, each=nrow(m))
	)
}

### Order
order.df = function(x, decreasing=TRUE, lim=c(1E+5, 10)) {
	if( ! is.null(lim)) {
		if(nrow(x) > lim[1] || ncol(x) > lim[2]) {
			msg = "Data frame is too big! Set lim = NULL to sort big data frames."
			stop(msg);
		}
	}
	order.s = function(...) order(..., decreasing=decreasing);
	id = do.call(order.s, x);
	return(id);
}

### Duplicates
countDuplicates = function(m, onlyDuplicates=FALSE) {
	count = rep(0, nrow(m));
	xd = aggregate(count ~ ., m, length);
	if(onlyDuplicates) {
		xd = xd[xd$count > 1, , drop=FALSE];
	}
	return(xd)
}

### Cut
# - cut by the Population Median;
# - compute proportions less/greater than population Median
#   in the various groups;
# - e = formula of type: lhs ~ groups;
# - N = total N;
cut.formula = function(e, data, FUN = median) {
	lhs = e[[2]];
	Mx_tmp = FUN(data[ , as.character(lhs)]);
	e[[2]] = as.call(list(as.symbol("<"), lhs, Mx_tmp));
	# e[[2]] = str2lang(paste0("(", lhs, " < ", Mx_tmp, ")"));
	FUNP = function(x) { s = sum(x) / length(x); c(s, 1 - s, length(x)); }
	dX.tbl = aggregate(e, data, FUNP)
	lvl = c("< Med", "> Med"); # c("LesserMed", "GreaterMed")
	len = ncol(dX.tbl);
	dX1 = dX.tbl; dX1$Type = factor(lvl[1], levels=lvl); dX1$Freq = dX.tbl[,len][,1]; dX1$N = dX.tbl[,len][,3];
	dX2 = dX.tbl; dX2$Type = factor(lvl[2], levels=lvl); dX2$Freq = dX.tbl[,len][,2]; dX2$N = dX.tbl[,len][,3];
	dX.tbl = rbind(dX1, dX2);
	dX.tbl = dX.tbl[, -len];
	dX.tbl
}


##################
##################

##################
### Formatting ###
##################

### Helper

# Argument matching
match.halign = function(justify, msg="Option for justify NOT supported!") {
	if(is.character(justify)) {
		id = pmatch(justify, c("left", "right", "center"));
		if(is.na(id) && justify == "centre") id = 3;
	}
	if(is.na(id)) stop(msg);
	return(id);
}
# String Operations
space.builder = function(nch, each=1, ch=" ") {
	chf = function(nch, each) rep(paste0(rep(ch, nch), collapse=""), each=each);
	sapply(nch, chf, each=each);
}
nchar.list = function(l) {
	lapply(l, nchar);
}
nchar.m = function(m) {
	nChD  = nchar(m);
	nDMax = apply(nChD, 2, max);
	return(list(max=nDMax, n=nChD));
}
pad.list = function(l, n, min=0, justify="right", ch=" ") {
	justify = match.halign(justify);
	nch = nchar.list(l);
	nsp = sapply(nch, function(n) max(n));
	nmx = pmax(nsp, min);
	ch0 = lapply(seq_along(nch), function(id)
		space.builder(nmx[[id]] - nch[[id]], each=1, ch=ch))
	pad.f = if(justify == 2) function(id) {
			paste0(ch0[[id]], l[[id]])
		} else if(justify == 1) function(id) {
			paste0(l[[id]], ch0[[id]])
		} else function(id) {
			pad.justify(l[[id]], nmx[[id]], nch[[id]], ch=ch);
		}
	l = lapply(seq_along(l), pad.f);
	attr(l, "nchar") = nmx;
	return(l);
}
pad.all = function(s, w, nch, justify="right", ch=" ") {
	justify = match.halign(justify);
	if(is.matrix(s)) w = rep(w, each=nrow(s));
	ch0 = sapply(seq_along(nch), function(id)
		space.builder(w[[id]] - nch[[id]], each=1, ch=ch));
	pad.f = if(justify == 2) function(id) {
			paste0(ch0[[id]], s[[id]])
		} else if(justify == 1) function(id) {
			paste0(s[[id]], ch0[[id]])
		} else function(id) {
			pad.justify(s[[id]], w[[id]], nch[[id]], ch=ch);
		}
	l = sapply(seq_along(s), pad.f);
	l.dim = dim(s);
	if( ! is.null(l.dim)) dim(l) = l.dim;
	attr(l, "nchar") = w;
	return(l);
}
pad.justify = function(s, nmax, nch, ch=" ") {
	if(missing(nch)) stop("nch: Not yet implemented!")
	nSpaces = nmax - nch;
	nLeft = nSpaces %/% 2; nRight = nSpaces - nLeft;
	mnCh = c(nLeft, nRight);
	ch0 = space.builder(mnCh, each=1, ch=ch);
	ch0 = matrix(ch0, ncol=2);
	paste0(ch0[,1], s, ch0[,2]);
}

# Merge 2 string matrices;
# Proper name: merge vs cbind?
merge.align = function(m1, m2, pos="Top", add.space=FALSE) {
	nr1 = nrow(m1); nr2 = nrow(m2);
	# TODO: "middle"-variants
	pos = if(is.numeric(pos)) pos else pmatch(pos, c("Top", "Bottom", "MiddleTop", "MiddleBottom"));
	# nchar
	getChars = function(m) {
		nch = attr(m, "nchar");
		if(is.null(nch)) nch = apply(m, 2, function(s) max(nchar(s)));
		return(nch);
	}
	nch1 = getChars(m1); nch2 = getChars(m2);
	# align
	if(nr1 > nr2) {
		if(add.space) {
			# add space to each cell of m2
			ch0 = space.builder(nch2, each = nr1 - nr2);
		} else {
			ch0 = matrix("", nrow = nr1 - nr2, ncol = ncol(m2));
		}
		m2 = if(pos == 1) rbind(m2, ch0) else rbind(ch0, m2);
	} else if(nr1 < nr2) {
		# add space to new rows of m1
		ch0 = space.builder(nch1, each = nr2 - nr1);
		m1 = if(pos == 1) rbind(m1, ch0) else rbind(ch0, m1);
	}
	m1 = cbind(m1, m2);
	attr(m1, "nchar") = c(nch1, nch2);
	return(m1);
}
rbind.align = function(m1, m2, justify="right", between=NULL) {
	# m1
	if(is.matrix(m1)) {
		nCh1  = nchar.m(m1);
		nMax1 = nCh1$max; nCh1 = nCh1$n;
	} else {
		nCh1 = nchar(m1); nMax1 = nCh1;
	}
	# m2
	nCh2  = nchar.m(m2);
	nMax2 = nCh2$max; nCh2 = nCh2$n;
	#
	wAll = pmax(nMax1, nMax2);
	m1 = pad.all(m1, w=wAll, nCh1, justify=justify);
	m2 = pad.all(m2, w=wAll, nCh2, justify=justify);
	if(is.null(between)) {
		m = rbind(m1, m2);
	} else {
		if(is.numeric(between)) {
			mB = space.builder(wAll, each = between, ch=" ");
		} else stop("Not yet implemented: \"between\"!");
		m = rbind(m1, mB, m2);
	}
	attr(m, "nchar") = wAll;
	return(m);
}
# Split names and align
split.names = function(names, min=0, extend=0, justify="right", pos="Top", split.ch = "\n",
			blank.rm=FALSE, detailed=TRUE, perl=TRUE) {
	# TODO: "Middle"
	justify = if(is.null(justify)) 1 else match.halign(justify);
	pos = if(is.null(pos)) 1 else pmatch(pos, c("Top", "Bottom", "MiddleTop", "MiddleBottom"));
	# Split strings
	str = strsplit(names, split.ch, perl=perl);
	if(blank.rm) str = lapply(str, function(s) s[nchar(s) > 0]);
	# nRows
	nr  = max(sapply(str, function(s) length(s)));
	# Width of each Column
	nch = sapply(str, function(s) max(nchar(s)));
	nch = pmax(nch, min);
	# Result
	ch0 = space.builder(nch, each=nr);
	mx  = matrix(ch0, nrow=nr, ncol=length(names));
	for(nc in seq(length(names))) {
		nrx = length(str[[nc]]); # current number of rows
		# Justifying
		nch.v = nchar(str[[nc]]);
		s = sapply(seq(nrx), function(nr) paste0(rep(" ", nch[[nc]] - nch.v[nr]), collapse=""));
		s = if(justify == 2) paste0(s, str[[nc]]) else if(justify == 1) paste0(str[[nc]], s)
			else {
				pad.justify(str[[nc]], nch[[nc]], nch.v, ch=" ");
			}
		if(pos == 1) {
			mx[seq(1, nrx), nc] = s;
		} else if(pos == 2) {
			mx[seq(nr + 1 - nrx, nr) , nc] = s;
		} else {
			# TODO: Middle-variants;
		}
	}
	if(detailed) attr(mx, "nchar") = nch;
	# Extend matrix: if option to extend;
	if(is.matrix(extend)) {
		mx = merge.align(mx, extend, pos=pos, add.space=TRUE);
	} else if(length(extend) > 1) {
		m.ext = matrix(space.builder(extend, each=nr), nrow=nr, ncol=length(extend));
		mx = cbind(mx, m.ext);
		attr(mx, "nchar") = c(nch, extend);
	} else if(extend > 0) {
		mx = cbind(mx, matrix("", nr=nr, ncol=extend));
		attr(mx, "nchar") = nch;
	}
	return(mx);
}
expand.labels = function(lst, default=" ", quote=FALSE) {
	len = lengths(lst);
	cpLensU = c(1, cumprod(len));
	cpLensD = rev(c(1, cumprod(rev(len))));
	y = NULL
	for (i in rev(seq_along(lst))) {
	    id = 1 + seq.int(from = 0, to = len[i] - 1) * cpLensD[i + 1L]
		ch0 = if(length(default) == 1) default else default[i];
	    tmp = rep(ch0, cpLensD[i])
	    tmp[id] = if(quote) charQuote(lst[[i]]) else lst[[i]];
	    y <- cbind(rep(tmp, times = cpLensU[i]), y)
	}
	y
}

### ftable with name splitting
# - this code should ideally replace format.ftable;
ftable2 = function(ftbl, print=TRUE, quote=FALSE, sep="|",
		justify="right", pos="Top", split.ch="\n", print.zero="-",
		method="auto", extend=TRUE, ...) {
	# Justify: the components of the argument
	if(is.list(justify)) {
		len = length(justify);
		if(len == 1) {
			justify = justify[[1]]; justify.lvl = justify; justify.num = "right";
		} else if(len == 2) {
			justify.lvl=justify[[2]]; justify = justify[[1]]; justify.num="right";
		} else {
			justify.lvl=justify[[2]]; justify.num = justify[[3]]; justify = justify[[1]];
		}
	} else {
		justify.lvl=justify; justify.num="right";
	}
	rvars = attr(ftbl, "row.vars");
	row.vars = names(rvars);
	cvars = attr(ftbl, "col.vars");
	col.vars = names(cvars);
	nV  = length(row.vars);
	ncv = length(cvars); # TODO
	# max width for each factor (all levels per factor);
	wL = sapply(nchar.list(rvars), max);
	nms = split.names(row.vars, min=wL, justify=justify, pos=pos,
		extend = 0, split.ch=split.ch);
	wL = attr(nms, "nchar")[seq(nV)];
	lvl = pad.list(rvars, min=wL, justify=justify.lvl);
	### Labels:
	ch0 = space.builder(wL, ch=" ");
	LBL = expand.labels(lvl, default = ch0, quote=quote);
	rownames(LBL) = NULL;
	
	### Data:
	DATA = sapply(ftbl, format, print.zero = print.zero, ...);
	dim(DATA) = dim(ftbl);
	# Column Vars
	ncc = ncv + sum(sapply(cvars, function(l) length(l))); # not yet used
	
	### Build ftbl
	LBL = rbind(nms, LBL);
	H0 = NULL;
	nrH = nrow(nms); between = NULL;
	method = pmatch(method, c("auto", "non.compact", "row.compact", "col.compact"));
	if(method == 3) {
		# "row.compact"
		if(nrH > 1) between = nrH - 1;
		nCh = nchar(col.vars);
		colVars = rbind(col.vars, space.builder(nCh, each = nrow(LBL) - 1));
		LBL = cbind(LBL, colVars);
	} else if(method == 4) {
		# "col.compact"
		# TODO
	} else if(method == 2) {
		between = nrH; # "non.compact"
		nH1 = sum(wL, (ncol(nms) - 1) * nchar(sep) );
		nCh = nchar(col.vars);
		colVars = space.builder(nCh, each = nrow(LBL));
		LBL = cbind(LBL, colVars);
		H0 = c(space.builder(nH1), col.vars);
	} else if(method == 1) {
		# "auto"
		H0 = space.builder(sum(wL, (ncol(nms) - 1) * nchar(sep) ));
		H0 = c(H0, col.vars); # TODO: all variables
		if(nrH > 1) between = nrH - 1;
	}
	
	### DATA
	DATA = rbind.align(unlist(cvars), DATA, between=between);
	if(method == 2) { H0 = c(H0, DATA[1,]); DATA = DATA[-1, ]; }
	
	### Build ftbl
	ftbl2 = cbind(LBL, DATA);
	if(print) {
		if( ! is.null(H0)) cat(H0, sep=c(rep(sep, length(H0) - 1), "\n"));
		cat(t(ftbl2), sep = c(rep(sep, ncol(ftbl2) - 1), "\n"))
	}
	invisible(ftbl2);
}

ftable2(ftbl, sep=" | ", zero.print="-", j="l", pos="Top", split.ch="\n", me="auto")


###############

### Encrypt IDs
encrypt = function(x, offset=0, isRandom=TRUE, DEBUG=TRUE) {
	# TODO: multiple columns in df;
	old.id = unique(x)
	len = length(old.id)
	if(DEBUG) print(len)
	
	new.id = seq(1+offset, len+offset)
	if(isRandom) {
		new.ids = sample(new.id, len)
	} else {
		new.ids = new.id
	}
	new.id = new.ids[match(x, old.id)]
	return(new.id)
}


##########################
##########################

### Formulas / Expressions
# - the Section has been moved to a separate file:
#   Tools.Formulas.R

