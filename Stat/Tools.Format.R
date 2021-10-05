########################
###
### Leonard Mada
### [the one and only]
###
### Format Tools
###
### draft v.0.1a


### Tools to Format Output


# shortcut
# source("Tools.Format.R")


###############
### History ###
###############


### draft v.0.1a:
# - moved formatting functions
#   from file: Tools.Data.R;


###########################
###########################

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

