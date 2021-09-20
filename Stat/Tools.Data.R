########################
###
### Leonard Mada
### [the one and only]
###
### Data Tools
###
### draft v.0.1k-fix2


### Tools to Process/Transform Data


###############
### History ###
###############


### draft v.0.1k - v.0.1k-fix2:
# - [fix] align = center;
# - [fix] proper argument: split.ch; [v.0.1k-fix]
# - [fix] added missing argument: pos; [v.0.1k-fix2]
### draft v.0.1j:
# - the Section on Formulas/Expressions
#   has been moved to a separate file:
#   Tools.Formulas.R;


#######################
#######################

### Groups / Aggregates

# TODO

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

# Helper
match.halign = function(justify, msg="Option for justify NOT supported!") {
	if(is.character(justify)) {
		id = pmatch(justify, c("right", "left", "center"));
		if(is.na(id) && justify == "centre") id = 3;
	}
	if(is.na(id)) stop(msg);
	return(id);
}
space.builder = function(nch, each=1, ch=" ") {
	chf = function(nch, each) rep(paste0(rep(ch, nch), collapse=""), each=each);
	sapply(nch, chf, each=each);
}
nchar.list = function(l) {
	lapply(l, nchar);
}
pad.list = function(l, n, min=0, justify="right", ch=" ") {
	justify = match.halign(justify);
	nch = nchar.list(l);
	nsp = sapply(nch, function(n) max(n));
	nmx = pmax(nsp, min);
	ch0 = lapply(seq_along(nch), function(id)
		space.builder(nmx[[id]] - nch[[id]], each=1, ch=ch))
	pad.f = if(justify == 1) function(id) {
			paste0(ch0[[id]], l[[id]])
		} else if(justify == 2) function(id) {
			paste0(l[[id]], ch0[[id]])
		} else function(id) {
			nSpaces = nmx[[id]] - nch[[id]];
			nLeft = nSpaces %/% 2; nRight = nSpaces - nLeft;
			mnCh = c(nLeft, nRight);
			ch0 = space.builder(mnCh, each=1, ch=ch);
			ch0 = matrix(ch0, ncol=2);
			paste0(ch0[,1], l[[id]], ch0[,2]);
		}
	l = lapply(seq_along(l), pad.f);
	attr(l, "nchar") = nmx;
	return(l);
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
		s = sapply(seq(nrx), function(nr) paste0(rep(" ", nch[[nc]] - nchar(str[[nc]][nr])), collapse=""));
		s = if(justify == 2) paste0(s, str[[nc]]) else paste0(str[[nc]], s);
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

### ftable with name splitting
# - this code should be ideally inside format.ftable;
ftable2 = function(ftbl, print=TRUE, quote=FALSE, sep="|",
		justify="right", justify.lvl=justify, pos="Top", extend=TRUE, split.ch="\n", ...) {
	rvars = attr(ftbl, "row.vars");
	row.vars = names(rvars);
	cvars = attr(ftbl, "col.vars");
	col.vars = names(cvars);
	nr  = length(row.vars);
	# Col columns
	ncc = length(col.vars) + sum(sapply(cvars, function(l) length(l)));
	nch = nchar(unlist(lapply(seq_along(cvars), function(id) c(col.vars[[id]], cvars[[id]]))));
	# max width for each factor (all levels per factor);
	w = sapply(nchar.list(rvars), max);
	nms = split.names(row.vars, min=w, justify=justify, pos=pos,
		extend = if(extend) nch else ncc, split.ch=split.ch);
	lvl = pad.list(rvars, min=attr(nms, "nchar")[seq(nr)], justify=justify.lvl);
	### format.ftbl
	# HACK: code should be ideally inside format.ftable!
	# - update width of factor labels;
	# - new width available in attr(nms, "nchar");
	tmp.lvl = lvl;
	names(tmp.lvl) = nms[1, seq_along(rvars)];
	attr(ftbl, "row.vars") = tmp.lvl; # use part of the name
	ftbl2 = format(ftbl, quote=quote, justify=justify, ...);
	# hack: insert the full names;
	ftbl2 = rbind(ftbl2[1,], nms, ftbl2[-c(1,2),]);
	if(print) {
		cat(t(ftbl2), sep = c(rep(sep, ncol(ftbl2) - 1), "\n"))
	}
	invisible(ftbl2);
}

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

