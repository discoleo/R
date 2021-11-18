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

order.df = function(x, decreasing=TRUE) {
	order.s = function(...) order(..., decreasing=decreasing);
	id = do.call(order.s, x);
	return(id);
}
### TODO: update use of sort.pm() everywhere!
sort.pm = function(p, xn=NULL, sort.coeff, xn2=NULL) {
	### Special Cols:
	# TODO: different approach;
	# - over xn: c(1,2,3,4) = Sum, Max, Min, MinNZ; (IF length(xn) > 1)
	# - over all: c(5,6,7,8,9) = SumAll, MaxAll, MinAll, MinNZAll, Coeff;
	isM = ( ! is.null(xn) && length(xn) > 1); # isMultiple
	if( ! is.null(xn2)) xn = c(xn, xn2);
	if(missing(sort.coeff)) {
		sort.coeff = if(isM) c(1,2, seq(10, length.out=length(xn)), 5,6)
			else if(is.null(xn)) c(1,2) else c(seq(6, length.out=length(xn)), 1,2);
	}
	pP = p[, - which(names(p) == "coeff"), drop=FALSE];
	summary.sort = function(p, FUN=sum) sapply(seq(nrow(p)), function(id) FUN(unlist(p[id, , drop=TRUE])));
	to.df = function(i, p, FUN) if(any(sort.coeff == i)) summary.sort(p, FUN) else rep(0, nrow(p));
	if(isM) {
		pP.pp = pP[, xn, drop=FALSE];
		s.df = to.df(1, pP.pp, sum);
		s.df = cbind(s.df, to.df(2, pP.pp, max));
		s.df = cbind(s.df, to.df(3, pP.pp, min));
		s.df = cbind(s.df, to.df(4, pP.pp, function(x) min(x[x != 0])));
		id0 = 4;
		s.df = cbind(s.df, to.df(1 + id0, pP, sum));
	} else {
		id0 = 0;
		s.df = to.df(1, pP, sum);
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
# used by: as.character.pm()
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
	p = cbind(p[,-leading, drop=FALSE], p[,leading, drop=FALSE]);
	return(p)
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
print.pm = function(...) {
	ch = as.character.pm(...);
	print(ch);
	invisible(ch);
}
print.p = function(...) {
	ch = as.character.pm(...);
	print(ch);
	invisible(ch);
}

as.character.pm = function(p, leading=NA, do.sort=TRUE, do.rev=FALSE, sort.order=TRUE,
		simplify.complex=TRUE, brackets.complex=TRUE) {
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
			p = cbind(p[,-leading, drop=FALSE], p[,leading, drop=FALSE]);
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
	isPlus = (Re(coeff) > 0) | (Re(coeff) == 0 & Im(coeff) >= 0);
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

toCoeff = function(p, x="x") {
	idx = match(x, names(p));
	if(idx < 0) stop(paste0("No variable ", x));
	px = p[,x]; p = p[, - idx, drop=FALSE];
	str = tapply(seq(nrow(p)), px, function(nr) print.p(p[nr,, drop=FALSE], leading=NA))
	str[nchar(str) == 0] = "1";
	# missing powers
	x.all = seq(0, max(px));
	p.all = rep("0", length(x.all));
	p.all[1 + sort(unique(px))] = str;
	return(p.all)
}
# Evaluate the coefficients using "..."
evalCoeff = function(p, xn="x", ...) {
	idx = match(xn, names(p));
	if(idx < 0) stop(paste0("No variable ", xn));
	px = p[,xn]; p = p[, - idx, drop=FALSE];
	if(ncol(p) > 1) {
		coeff = tapply(seq(nrow(p)), px, function(nr) eval.pm(p[nr,, drop=FALSE], ...));
	} else coeff = p$coeff;
	# missing powers
	x.all = seq(0, max(px));
	p.all = rep(0, length(x.all));
	p.all[1 + sort(unique(px))] = coeff;
	p.all = rev(p.all);
	return(p.all);
}
coef.pm = function(p, xn="x", descending=TRUE) {
	if(ncol(p) > 2) warning("Multi-variable polynomial!");
	p = aggregate0.pm(p[, c(xn, "coeff"), drop=FALSE]);
	p = reduce.pm(p);
	# missing powers
	p.all = rep(0, max(p[, xn]) + 1);
	p.all[p[, xn] + 1] = p$coeff;
	if(descending) p.all = rev(p.all);
	return(p.all);
}
### Print
print.coeff = function(p, x="x") {
	p = rev(toCoeff(p, x));
	last = tail(p, 1);
	sapply(head(p, -1), function(p) cat(paste(p, ",\n", sep="")));
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

