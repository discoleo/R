########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Format & Print


### fast load:
# - is automatically loaded in: Polynomials.Helper.R;
# source("Polynomials.Helper.Format.R")


######################

### General Functions

# TODO: analyse efficiency;
cut.character.int = function(n, w, extend=2) {
	if(w == 0) {
		len = length(n);
		pos = rbind(seq(len), seq(len));
		return(pos);
	}
	n = n + extend;
	cumlen = cumsum(n);
	max = tail(cumlen, 1) %/% w + 1;
	pos = cut(cumlen, seq(0, max) * w);
	count = rle(as.numeric(pos))$lengths;
	pos = cumsum(count);
	posS = pos[ - length(pos)] + 1;
	posS = c(1, posS);
	pos = rbind(posS, pos);
	return(pos);
}
cut.character.int.old = function(n, w) {
	ncm = cumsum(n);
	nwd = ncm %/% w;
	count = rle(nwd)$lengths;
	pos = cumsum(count);
	posS = pos[ - length(pos)] + 1;
	posS = c(1, posS);
	pos = rbind(posS, pos);
	return(pos);
}

# Basic function to order a df;
order.df = function(x, decreasing=TRUE) {
	order.s = function(...) order(..., decreasing=decreasing);
	id = do.call(order.s, x);
	return(id);
}

### Specific Functions

### TODO: update use of sort.pm() everywhere!
# TODO: do.max;
# do.sum = 0: NO sum; MAX has priority, then individual powers;
#    x^4 > x^3*y^4 > y^4 > x^3*y*z > x^3*z^3;
# do.sum = 1: sum has priority, then MAX and individual powers;
#    x^4 > y^4 > x^3*y > x^2*y^2 > x^3;
# do.sum = 2: MAX has priority, then sum and individual powers;
#    x^4 > x^3*y*z > x^3*z^2 > x^3*y; (only 1st MAX counts)
sort.pm = function(p, xn=NULL, xn2=NULL, do.sum=1, sort.coeff) {
	### Special Cols:
	# TODO: different approach;
	# - over xn: c(1,2,3,4) = Sum, Max, Min, MinNZ; (IF length(xn) > 1)
	# - over all: c(5,6,7,8,9) = SumAll, MaxAll, MinAll, MinNZAll, Coeff;
	len = if(is.null(xn)) 0 else length(xn);
	isM = (len > 1); # isMultiple
	xnM = xn; # Order: (Sum, Max) => x^2, y^2, x*y;
	if( ! is.null(xn2)) xn = c(xn, xn2);
	if(missing(sort.coeff)) {
		idSort =
			if(do.sum == 1) { c(1,2); }
			else if(do.sum == 0) c(2) else c(2,1);
		# Sort priorities:
		sort.coeff = if(isM) c(idSort, seq(10, length.out=length(xn)), 5,6)
			else if(is.null(xn)) idSort
			else c(seq(6, length.out=length(xn)), idSort);
	}
	# Check if Polynomial:
	idCoeff = which(names(p) == "coeff");
	if(length(idCoeff) != 1) stop("Missing Coefficients!");
	pP = p[, - idCoeff, drop=FALSE];
	#
	summary.sort = function(p, FUN=sum) sapply(seq(nrow(p)), function(id) FUN(unlist(p[id, , drop=TRUE], recursive=FALSE)));
	to.df = function(i, p, FUN) if(any(sort.coeff == i)) summary.sort(p, FUN) else rep(0, nrow(p));
	if(isM) {
		# sum only over xnM:
		pP.pp = pP[, xnM, drop=FALSE];
		s.df = to.df(1, pP.pp, sum);
		s.df = cbind(s.df, to.df(2, pP.pp, max));
		s.df = cbind(s.df, to.df(3, pP.pp, min));
		s.df = cbind(s.df, to.df(4, pP.pp, function(x) min(x[x != 0])));
		id0 = 4;
		# sum over all variables:
		s.df = cbind(s.df, to.df(1 + id0, pP, sum));
	} else {
		id0 = 0;
		s.df = to.df(1, pP, sum);
		s.df = as.data.frame(s.df);
	}
	s.df = cbind(s.df, to.df(2 + id0, pP, max));
	s.df = cbind(s.df, to.df(3 + id0, pP, min));
	s.df = cbind(s.df, to.df(4 + id0, pP, function(x) min(x[x != 0]) ));
	s.df = cbind(s.df, to.df(5 + id0, p[, "coeff", drop=FALSE], function(x) - abs(x) ));
	if( ! is.null(xn)) s.df = cbind(s.df, pP[, xn, drop=FALSE]);
	s.df = s.df[, sort.coeff, drop=FALSE];
	id = order.df(s.df, decreasing=TRUE);
	p = p[id,]; rownames(p) = seq(nrow(p));
	return(p)
}
# set Coeff Column as last column;
# TODO: more;
sortColumns.pm = function(p) {
	idc = match("coeff", names(p));
	p = cbind(p[ , - idc, drop=F], coeff=p[, idc]);
	return(p)
}
### Order Variables
# - sorts variables in specified order;
# - returns the actual polynomial;
orderVars.pm = function(p, xn, last=TRUE, warn=TRUE) {
	stop("Deprecated! Use: orderVars (non-generic).")
}
orderVars = function(p, xn, last=TRUE, warn=TRUE, coeff.last=TRUE) {
	sort.pm.vars(p, xn=xn, last=last, warn=warn, coeff.last=coeff.last);
}
sort.pm.vars = function(p, xn, last=TRUE, warn=TRUE, coeff.last=TRUE) {
	if(last && coeff.last) {
		hasCoeff = match("coeff", xn);
		if(is.na(hasCoeff)) {
			xn = c(xn, "coeff");
		}
	}
	ids = match(xn, names(p));
	isNA = is.na(ids);
	if(any(isNA)) {
		if(warn) warning("Variables: ", xn[isNA], " are NOT present!");
		xn = xn[ ! isNA];
	}
	id0 = match(names(p)[ids], xn);
	#
	if(last) {
		p = cbind(p[ , -ids, drop=FALSE], p[ , ids[id0], drop=FALSE]);
	} else {
		p = cbind(p[ , ids[id0], drop=FALSE], p[ , -ids, drop=FALSE]);
	}
	return(p);
}

### proper Order of variables
sort.pm.proper = function(p, xn = c("b", "R", "E", "S"), warn=TRUE, do.grep=TRUE) {
	xnr = if(do.grep) paste0("^", xn) else xn;
	nms = names(p);
	# IDs of variables
	ids = integer(0);
	for(nm in xnr) {
		# TODO: do.grep = FALSE;
		id = grep(nm, nms);
		if(length(id) == 0) {
			if(warn) warning("Variable: ", nm, " not found!");
			next;
		}
		id  = id[order(nms[id])];
		ids = c(ids, id);
	}
	idc = match("coeff", nms);
	if( ! idc %in% ids) ids = c(ids, idc);
	if(length(ids) == ncol(p)) {
		p = p[, ids, drop=FALSE];
		return(p);
	}
	tmp = p[, - ids, drop=FALSE];
	p = cbind(tmp, p[, ids, drop=FALSE]);
	return(p)
}

# used by: as.character.pm()
# sort.order = default decreasing;
# do.rev = reverse order of multiple leading variables:
#   (x,y,z) => (z,y,x);
sort.simple.pm = function(p, leading=1, do.rev=FALSE, sort.order=TRUE) {
	if(length(leading) == 1) {
		p = p[order(p[, leading], decreasing=sort.order), , drop=FALSE];
	} else {
		coeff.df = data.frame(vs=apply(p[, leading, drop=FALSE], 1, sum));
		coeff.df = cbind(coeff.df, (p[, leading]));
		order.s = function(...) order(..., decreasing=sort.order);
		id = do.call(order.s, coeff.df);
		p = p[id,];
		if(do.rev) leading = rev(leading);
	}
	orderVars = ! is.na(do.rev);
	if(orderVars) {
		p = cbind(p[ , -leading, drop=FALSE], p[ , leading, drop=FALSE]);
	}
	return(p)
}

as.pm.lead = function(p, xn, warn=TRUE) {
	# is moved as last column;
	id = match(xn, names(p));
	if(is.na(id)) {
		if(warn) warning("Variable not found!");
		return(p);
	}
	p = cbind(p[, -id, drop=FALSE], p[ , id, drop=FALSE]);
	return(as.pm(p));
}
as.pm.first = function(p, xn, warn=TRUE) {
	# is moved as first column;
	id = match(xn, names(p));
	if(is.na(id)) {
		if(warn) warning("Variable not found!");
		return(p);
	}
	p = cbind(p[ , id, drop=FALSE], p[, -id, drop=FALSE]);
	return(as.pm(p));
}

#############

#############
### Print ###
#############

# Print multi-variable Poly
print.monome = function(name, p) {
	v = p[,name];
	v.r = rep("", length(v));
	isCoeff = v != 1 & v != 0;
	v.r[isCoeff] = paste0(name, "^", v[isCoeff]);
	v.r[v == 1] = name;
	return(v.r);
}
print.pm = function(..., print=TRUE) {
	ch = as.character.pm(...);
	if(print) { cat(ch); cat("\n"); }
	invisible(ch);
}
print.p = function(..., print=TRUE) {
	ch = as.character.pm(...);
	if(print) { cat(ch); cat("\n"); }
	invisible(ch);
}
# print as data.frame;
print.df = function(p, n=100) UseMethod("print.df");
print.df.pm = function(p, n=100) {
	if(is.null(n) || is.na(n)) {
		print.data.frame(p);
		return(invisible());
	}
	head(as.data.frame(p), n=n);
}

# - leading = "leading" variable, printed at the end of the monomials;
# - do.sort = sort monomials based on "leading" variables;
# - do.rev = order "leading" variables in reverse order;
# - sort.order = default descending;
as.character.pm = function(p, leading=NA, do.sort=TRUE, do.rev=FALSE, sort.order=TRUE,
		simplify.complex=TRUE, brackets.complex=TRUE) {
	if(inherits(p, "pm.div")) {
		if(nrow(p$Rem) > 0) warning("The Remainder of division is NOT printed!");
		p = p$Rez;
	}
	if(nrow(p) == 0) return("");
	### Var order
	isNA = all(is.na(leading));
	if( ! isNA && ! is.numeric(leading)) leading = match(leading, names(p));
	if(any(is.na(leading))) {
		if( ! isNA) warning("Sort var does NOT exist!");
		leading = leading[ ! is.na(leading)];
	}
	if(length(leading) > 0) {
		if(do.sort) {
			p = sort.simple.pm(p, leading=leading, do.rev=do.rev);
		} else {
			if(do.rev) leading = rev(leading);
			p = cbind(p[ , -leading, drop=FALSE], p[ , leading, drop=FALSE]);
		}
	}
	###
	id.coeff = match("coeff", colnames(p));
	coeff = p[,id.coeff]; p = p[, - id.coeff, drop=FALSE];
	p.str = sapply(colnames(p), print.monome, p=p);
	# print(p.str)
	paste.nonempty = function(str, collapse="*") {
		str = str[nchar(str) > 0]
		paste(str, collapse=collapse)
	}
	if( ! is.null(dim(p.str))) p.str = apply(p.str, 1, paste.nonempty)
	else p.str = paste.nonempty(p.str);
	# Sign: 0 treated as "+"
	isPlus = if(inherits(coeff, c("bigz", "bigq"))) (coeff > 0)
		else (Re(coeff) > 0) | (Re(coeff) == 0 & Im(coeff) >= 0);
	sign.str = ifelse(isPlus, " + ", " - ");
	sign.str[1] = if(isPlus[1]) "" else "- ";
	# Complex numbers
	coeffPlus = as.abs.complex(coeff, rm.zero = simplify.complex, coupled = brackets.complex);
	coeff.str = format.complex.pm(coeffPlus, rm.zero = simplify.complex, brackets = brackets.complex);
	# Coeff == 1
	hasCoeff = (coeffPlus != 1 & nchar(p.str) > 0); # TODO: verify if fixed!
	isB0 = (nchar(p.str) == 0);
	p.str[hasCoeff] = paste(coeff.str[hasCoeff], p.str[hasCoeff], sep = "*");
	p.str[isB0] = coeff.str[isB0]; # [should be fixed]: ERROR "+ b0*";
	return(paste(sign.str, p.str, sep="", collapse=""));
}
# coupled: "-2 + 3i" => - "(2 - 3i)";
# de-coupled: "-2 + 3i" => - "2 + 3i";
as.abs.complex = function(x, rm.zero=TRUE, coupled=TRUE) {
	if(inherits(x, c("bigz", "bigq"))) return(abs(x));
	isNegativ = Re(x) < 0;
	if(is.numeric(x)) {
		x[isNegativ] = - x[isNegativ];
		return(x);
	}
	# Complex
	x[isNegativ] = if(coupled) { - x[isNegativ]; }
		else complex(re = - Re(x[isNegativ]), im = Im(x[isNegativ]));
	if(rm.zero) {
		isImNegativ = (Re(x) == 0 & Im(x) < 0);
		x[isImNegativ] = complex(0, im = - Im(x[isImNegativ]));
	}
	return(x);
}
# Note:
# format.complex: breaks format.default!
format.complex.pm = function(x, sign.invert=FALSE, rm.zero=TRUE, brackets=TRUE, i.ch="i") {
	if(sign.invert) {
		x = as.abs.complex(x, rm.zero = rm.zero, coupled = ! is.null(brackets));
	}
	coeff.str = as.character(x);
	if(inherits(x, c("bigz", "bigq"))) return(coeff.str);
	if(rm.zero) {
		isReZero = ((Re(x) == 0) & (Im(x) != 0));
		coeff.str[isReZero] = paste0(Im(x[isReZero]), i.ch);
	}
	if(brackets && is.complex(x)) {
		coeff.str = paste0("(", coeff.str, ")");
	}
	return(coeff.str);
}

### Other

### Convert to Coefficients
# - as list of polynomials or of numeric values;
# - the list is in descending order;
as.coeff.pm = function(p, xn) {
	idv = match(xn, names(p));
	# b0:
	if(is.na(idv)) {
		warning("Variable ", xn, " not found!");
		return(list(p));
	}
	#
	if(ncol(p) == 2) {
		coeff = coef.pm(p, xn=xn, descending=TRUE);
		return(as.list(coeff));
	}
	#
	pows = seq(max(p[ , idv]), 0, by=-1);
	rez  = lapply(pows, function(pow) {
		pb = p[p[, idv] == pow, - idv, drop=FALSE];
		if(nrow(pb) == 0) return(0);
		return(pb);
	})
	return(rez);
}

### Convert to Coefficients: as string;
# TODO: check everywhere that x is replaced with xn;
toCoeff = function(p, xn="x", decreasing=TRUE, print=TRUE, addComments=FALSE,
		sep=NULL, WIDTH=80) {
	idx = match(xn, names(p));
	if(idx < 0) stop(paste0("No variable ", xn));
	px = p[,xn]; p = p[, - idx, drop=FALSE];
	str = tapply(seq(nrow(p)), px, function(nr) as.character.pm(p[nr,, drop=FALSE], leading=NA))
	str[nchar(str) == 0] = "1";
	# missing powers
	x.all = seq(0, max(px));
	p.all = rep("0", length(x.all));
	p.all[1 + sort(unique(px))] = str;
	if(decreasing) p.all = rev(p.all);
	class(p.all) = c("pm.coeff", class(p.all));
	attr(p.all, "xn") = list(xn = xn, isDesc=decreasing);
	### Print:
	if(print) {
		cat.pm.coeff(p.all, sep=sep, w=WIDTH, addComments=addComments);
		# return invisibly: coefficients are already printed;
		return(invisible(p.all));
	}
	return(p.all)
}
cat.pm.coeff = function(p, sep=NULL, w=60, addComments=FALSE) {
	if(is.null(sep)) {
		LEN = length(p);
		if(LEN <= 1) { sep = "\n"; }
		else {
			pos = cut.character.int(nchar(p), w=w);
			nc  = ncol(pos);
			if(addComments) {
				xn = attr(p, "xn");
				if(is.null(xn)) { addComments = FALSE; }
				else xn = xn$xn; # TODO: isDesc;
			}
			#
			if(nc > 1) for(id in seq(nc - 1)) {
				if(addComments) {
					npow = pos[2, nc] - pos[1, id];
					if(npow > 1) { spow = "^"; }
					else { spow = ""; npow = ""; }
					cat(paste0("# ", xn, spow, npow, "\n"), sep="");
				}
				slen = pos[2, id] - pos[1, id];
				# "" vs ",\n": BUG in cat() inside FOR loop?
				xsep = c(rep(", ", slen), "");
				cat(p[seq(pos[1, id], pos[2, id])], sep = xsep); cat(",\n");
			}
			if(addComments) cat("# B0", sep="\n");
			slen = pos[2, nc] - pos[1, nc];
			cat(p[seq(pos[1, nc], pos[2, nc])], sep = c(rep(", ", slen), "\n"));
			return();
		}
	}
	xsep = c(rep(sep, length(p) - 1), "\n");
	cat(p, sep = xsep);
}

# Evaluate the coefficients using "..."
evalCoeff = function(p, xn="x", ...) {
	idx = match(xn, names(p));
	if(is.na(idx)) stop(paste0("No variable ", xn));
	px = p[,xn]; p = p[, - idx, drop=FALSE];
	if(ncol(p) > 1) {
		coeff = tapply(seq(nrow(p)), px, function(nr) eval.pm(p[nr,, drop=FALSE], ...));
		px = sort(unique(px));
	} else {
		coeff = p$coeff;
		if(any(duplicated(px))) stop("TODO: Duplicated powers!");
	}
	# Missing powers
	x.all = seq(0, max(px));
	p.all = rep(0, length(x.all));
	p.all[1 + px] = coeff;
	# Ascending order:
	p.all = rev(p.all);
	return(p.all);
}
# returns only the numeric coefficients
coef.pm = function(p, xn=NULL, pow=NULL, descending=TRUE) {
	if(is.null(xn)) {
		if(ncol(p) > 2) stop("Missing variable name!");
		idc = match("coeff", names(p));
		if(is.na(idc)) stop("Missing coefficients!");
		xn = names(p)[ - idc];
	} else if(ncol(p) > 2) warning("Multi-variable polynomial!");
	p = aggregate0.pm(p[, c(xn, "coeff"), drop=FALSE]);
	p = reduce.pm(p);
	# Missing powers
	if(is.null(pow)) {
		p.all = rep(0, max(p[, xn]) + 1);
		p.all[p[, xn] + 1] = p$coeff;
		if(descending) p.all = rev(p.all);
	} else {
		p.all = rep(0, length(pow));
		id = match(pow, p[, xn]);
		isNotNA = ! is.na(id);
		id = id[isNotNA];
		p.all[isNotNA] = p$coeff[id];
		# Descending: ???
	}
	return(p.all);
}
### Print
print.coeff = function(p, x="x") {
	p = rev(toCoeff(p, x));
	last = tail(p, 1);
	err = sapply(head(p, -1), function(p) cat(paste(p, ",\n", sep="")));
	cat(paste(last, "\n", sep=""));
	invisible(p);
}
print.pcoeff = function(l, print=TRUE, strip=NULL, len=10) {
	nlast = length(l);
	lsep = rep(", ", nlast);
	lsep[nlast] = ""; # tail(lsep, 1) = ""; # DOES NOT function!
	lsep[seq(len+1, nlast-1, by=len)] = ",\n";
	l.str = paste0(l, lsep, collapse="");
	if( ! is.null(strip)) {
		l.str = gsub(paste0("[ *]*+", strip, "\\^[0-9]++"), "", l.str, perl=TRUE);
	}
	if(print) { cat(l.str); cat("\n"); }
	return(invisible(l.str));
}

