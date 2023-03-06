########################
###
### Leonard Mada
### [the one and only]
###
### Data Tools
###
### draft v.0.2d


### Tools to Process/Transform Data


### Note:
# - the Formatting helper functions were moved
#   to file Tools.Format.R;

### fast load:
# source("Tools.Data.R");

###############
### History ###
###############


### draft v.0.2c:
# - [new] Function: as.na;
### draft v.0.2b:
# - mid.factor: compute mid-points of an interval;
### draft v.0.2a:
# - moved section with Formatting functions
#   to new file: Tools.Format.R;
### draft v.0.1q - v.0.1q-ref:
# - recode factor levels: with fail-safe provisions;
# - reordered sections; [v.0.1q-reorder]
# - [refactor] argument name: more.warn; [v.0.1q-ref]
### draft v.0.1p:
# - basic implementation of a function to encode
#   numeric values as a factor: as.factor.df();
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


######################
######################

### DF Data Transforms

### NA
as.na = function(x, cols, code=99) {
	if(inherits(x, "data.frame")) {
		for(nm in cols) {
			x[x[, nm] == code, nm] = NA;
		}
		return(x);
	}
	stop("Not yet implemented!");
}

### Encode: Numeric => Factor
as.factor.df = function(x, vars, name = c("Lvl ", ""), ordered=TRUE) {
	if( ! inherits(x, "data.frame"))
		stop("Error: x must be a data frame!");
	idc = match(vars, names(x));
	if(any(is.na(idc))) stop(
		paste0("Error: column ", vars[is.na(idc)][1], " does not exist!"));
	if(is.null(dim(name))) {
		if(length(name) == 1) name = c(name, "");
	}
	for(id in idc) {
		x[, id] = factor(x[, id], ordered=ordered);
		# TODO: matrix with names;
		levels(x[, id]) = paste0(name[1], levels(x[, id]), name[2]);
	}
	return(x);
}

### Recode factors
recode.df.factor = function(x, var.name, ..., more.warn = "Stop") {
	x.f = x[ , var.name, drop=FALSE];
	if(ncol(x.f) != 1) stop("Error: can process only 1 column!");
	x.f = x.f[ , 1, drop=TRUE];
	if( ! is.factor(x.f)) stop("Error: can process only factors!")
	lvl.old = levels(x.f);
	### Levels
	e = substitute(c(...));
	len = length(e);
	ech = as.character(e)[-1]; # new levels
	lvl.nms.old = names(e)[-1]; # old levels
	len = len - 1; idDot = -1;
	hasDot = FALSE;
	for(id in seq(len)) {
		if(ech[[id]] == ".") {
			if(hasDot) stop("Error: duplicate dot!");
			idDot  = id;
			hasDot = TRUE;
		}
	}
	if(hasDot) {
		len = len - 1;
		# ">=" vs ">": What is best for c(same.levels, .) ?
		if(len >= length(lvl.old)) {
			more.warn = pmatch(more.warn, c("Stop", "Warn", "Ignore"));
			if(is.na(more.warn)) stop("Error: option for more.warn not supported!");
			if(more.warn == 1) stop("Error: too many levels!");
			if(more.warn == 2) warning("Too many levels!");
		}
		ech = ech[-idDot]; lvl.nms.old = lvl.nms.old[-idDot];
	} else if(len != length(lvl.old))
		stop("Error: mismatch in number of levels!");
	idl = match(lvl.nms.old, lvl.old);
	isNA = is.na(idl);
	# TODO: handle new levels
	# could use: . = "New Level";
	if(any(isNA)) {
		stop(paste0("Error: levels", paste0(lvl.nms.old[isNA], collapse=", "), " do NOT exist!"));
	}
	lvl.new = lvl.old;
	lvl.new[idl] = ech;
	levels(x.f) = lvl.new;
	x[ , var.name] = x.f;
	return(x);
}


### Factors

### Middle of an Interval
mid.factor = function(x, inf.to = NULL, split.str=",") {
	lvl0 = levels(x); lvl = lvl0;
	lvl = sub("^[(\\[]", "", lvl);
	lvl = sub("[])]$", "", lvl); # tricky Regex;
	lvl = strsplit(lvl, split.str);
	lvl = lapply(lvl, function(x) as.numeric(x));
	if( ! is.null(inf.to)) {
		FUN = function(x) {
			if(any(x == Inf)) 1
			else if(any(x == - Inf)) -1
			else 0;
		}
		whatInf = sapply(lvl, FUN);
		# TODO: more advanced;
		lvl[whatInf == -1] = inf.to[1];
		lvl[whatInf ==  1] = inf.to[2];
	}
	mid = sapply(lvl, mean);
	lvl = data.frame(lvl=lvl0, mid=mid);
	merge(data.frame(lvl=x), lvl, by="lvl");
}

### Grid / Hierarchical Data

# freq = counts for level 2;
expand.hierarchy = function(freq, labels=NULL) {
	n = length(freq);
	h1 = unlist(lapply(seq(n), function(i) rep(i, freq[i])));
	h2 = unlist(lapply(seq(n), function(i) if(freq[i] == 0) NULL else seq(freq[i])));
	r = data.frame(h1, h2);
	if( ! is.null(labels)) {
		colnames(r) = labels;
	}
	return(r);
}


##################
##################

### DF Transforms

### Rename
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

#######################

### Groups / Aggregates

# TODO

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


### Helper functions:
# - moved to file: Tools.Format.R;
source("Tools.Format.R")


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


if(FALSE) {
	# not run
	ftable2(ftbl, sep=" | ", zero.print="-", j="l", pos="Top", split.ch="\n", me="auto")
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

