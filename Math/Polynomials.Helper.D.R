########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Derivation
###
### draft v.0.1d


### fast load:
# source("Polynomials.Helper.D.R")

### requires:
# - but is already loaded inside;
#   source("Polynomials.Helper.R")

######################

###############
### History ###
###############


### draft v.0.1d:
# - basic Integration;
### draft v.0.1c - v.0.1c-fix:
# - split.pm: prepare polynomial fractions for Integration;
# - [fix] absent variable name;
### draft v.0.1b:
# - various new methods:
#   generate various types of ODEs, e.g.
#   dp.log.pm, etc;
### draft v.0.1a:
# - moved Section on Derivation
#   from Polynomials.Helper.R;


#####################
#####################

### General Methods

sort.dpm = function(p, y="y", x="x") {
	# TODO: "dny"
	nms = if(length(y) > 1) y else paste0(c("d2", "d", ""), y);
	nms = c(nms, x);
	# Valid names:
	idSort = nms %in% names(p);
	nms = nms[idSort];
	if(length(nms) == 0) warning("No valid names!");
	# hard-coded: start.id = 10
	p = sort.pm(p, nms, sort.coeff=seq(10, length.out=length(nms)));
	nms = c(names(p)[ ! names(p) %in% nms], rev(nms));
	p = p[, nms];
	return(p);
}
print.dpm = function(p, y="y", x="x", do.sort=TRUE) {
	if(do.sort) p = sort.dpm(p, y=y, x=x);
	print.pm(p, do.sort=FALSE, leading=NA);
}

### D

# D( p(x) )
dp.pm = function(p, by="x", xn=by) {
	nc = match(xn, names(p));
	if(is.na(nc)) return(0);
	p = p[(p[, nc] != 0), , drop=FALSE];
	if(nrow(p) == 0) return(0);
	p$coeff = p$coeff * p[, nc];
	p[, nc] = p[, nc] - 1;
	return(p);
}
dnp.pm = function(p, n=2, by="x", xn=by) {
	for(id in seq(n)) {
		p = dp.pm(p, xn=xn);
		if(is.numeric(p)) break;
	}
	return(p);
}
# Dn( p(x), n = max(p$x) );
dp.pm.all = function(p, by, reduce=TRUE, warn=TRUE) {
	if(missing(by)) stop("Missing variable!");
	if(length(by) > 1) {
		warning("More than 1 variable! Processing sequentially.");
		for(xn in by) {
			p = dp.pm.all(p, by=xn, reduce=reduce);
		}
		return(p);
	}
	# Fast technique:
	id = match(by, names(p));
	if(is.na(id)) {
		if(warn) warning("variable ", by, " not present!");
		return(data.frame(coeff=0));
	}
	xmax = max(p[, id]);
	if(xmax == 0) {
		if(warn) warning("variable ", by, " not present!");
		return(data.frame(coeff=0));
	}
	p = p[p[, id] == xmax, -id, drop=FALSE];
	p = drop.pm(p);
	# Reduce using the gcd:
	if(reduce) {
		div = gcd.vpm(p, xgcd=p$coeff[1]);
		if(div > 1) p$coeff = p$coeff / div;
	} else {
		if(inherits(p$coeff, c("bigz", "bigq"))) {
			f = factorialZ(xmax);
		} else f = factorial(xmax);
		p$coeff = p$coeff * f;
	}
	return(p);
}

### Specific Derivatives

# D( p$Poly * exp(p$Exp) )
# p = list(Poly, Exp);
dp.exp.pm = function(p, xn="x") {
	pr = mult.pm(dp.pm(p$Exp, xn), p$Poly);
	pr = sum.pm(pr, dp.pm(p$Poly, xn));
	p$Poly = pr;
	return(p);
}
# D( p$C * log(p$Log) )
dp.log.pm = function(p, xn="x") {
	pB0 = mult.pm(dp.pm(p$Log, xn), p$C);
	pPP = mult.pm(p$Log, dp.pm(p$C, xn));
	p$C = pPP; p$B0 = pB0; p$Div = p$Log;
	return(p);
}
# D( p1 / pdiv ])
# TODO: consistent output-names;
dp.div.pm = function(p1, pdiv, xn="x") {
	if(is.numeric(p1) || is.complex(p1)) {
		r = mult.pm(dp.pm(pdiv, xn=xn), - p1);
		div = mult.pm(pdiv, pdiv);
		return(list(p = r, div = div));
	}
	r  = mult.pm(pdiv, dp.pm(p1, xn=xn));
	r2 = mult.pm(p1, dp.pm(pdiv, xn=xn));
	r = diff.pm(r, r2);
	r = list(p=r, div = mult.pm(pdiv, pdiv));
	return(r)
}
# D( p1*sin(T(x)) + p2*cos(T(x)) )
dp.trig.pm = function(plst, pT, xn="x", trig.order="sin") {
	if(inherits(plst, "pm")) {
		p1 = plst;
		p2 = NULL;
	} else if(is.list(plst)) {
		isD = ! is.null(attr(plst, "D"));
		if(isD) {
			p1 = plst$C1; p2 = plst$C2;
			if( ! missing(pT)) warning("Using pT included in the list!")
			pT = plst$Trig;
			trig.order = attr(plst, "trig.order");
		} else {
			p  = extract2.pm(plst, "p list");
			p1 = p[[1]]; p2 = p[[2]];
		}
	}
	# Order: p1 * sin, then p2 * cos;
	trig.order.id = match(trig.order, c("sin", "cos"));
	if(is.null(trig.order.id)) stop("Invalid order: p[[1]]*sin(T) + p[[2]]*cos(T);");
	# Result
	dT = function(p, pT, type) {
		r = dp.pm(pT, xn=xn);
		if(type == 2) r$coeff = - r$coeff;
		r = mult.pm(r, p);
		return(r);
	}
	next.type = function(id) if(id == 1) 2 else 1;
	C1 = dp.pm(p1, xn=xn);
	C2 = dT(p1, pT, type=trig.order.id);
	if( ! is.null(p2)) {
		C1 = sum.pm(C1, dT(p2, pT, type = next.type(trig.order.id)));
		C2 = sum.pm(C2, dp.pm(p2, xn=xn));
	}
	rlst = list(C1=C1, C2=C2, Trig=pT);
	attr(rlst, "trig.order") = trig.order;
	attr(rlst, "D") = "Trig";
	return(rlst);
}

# D( p1*sin(T0(x) + log(T1(x))) + p2*cos(T0(x) + log(T1(x))) )
dp.trigLog.pm = function(plst, pT, xn="x", trig.order="sin") {
	if(inherits(plst, "pm")) {
		p1 = plst;
		p2 = NULL;
		pT = extract2.pm(pT, "pT list");
		pT1 = pT[[1]]; pT0 = pT[[2]];
	} else if(is.list(plst)) {
		isD = ! is.null(attr(plst, "D")); # TODO: type;
		if(isD) {
			p1 = plst$C1; p2 = plst$C2;
			if( ! missing(pT)) warning("Using pT included in the list!");
			pT1 = plst$TrigLog;
			pT0 = plst$Trig0;
			trig.order = attr(plst, "trig.order");
		} else {
			p  = extract2.pm(plst, "p list");
			p1 = p[[1]]; p2 = p[[2]];
			pT = extract2.pm(pT, "pT list");
			pT1 = pT[[1]]; pT0 = pT[[2]];
		}
	}
	# Order: p1 * sin, then p2 * cos;
	trig.order.id = match(trig.order, c("sin", "cos"));
	if(is.null(trig.order.id)) stop("Invalid order: p[[1]]*sin(T) + p[[2]]*cos(T);");
	# Result
	dT = function(p, type) {
		r = dp.pm(pT1, xn=xn);
		if( ! is.null(pT0)) r = sum.pm(r, mult.pm(pT1, dp.pm(pT0, xn=xn)));
		if(type == 2) r$coeff = - r$coeff;
		r = mult.pm(r, p);
		return(r);
	}
	next.type = function(id) if(id == 1) 2 else 1;
	C1 = mult.pm(pT1, dp.pm(p1, xn=xn)); # T1 * D(P1)
	C2 = dT(p1, type=trig.order.id); # P1*(T1*D(T0) + D(T1))
	if( ! is.null(p2)) {
		C1 = sum.pm(C1, dT(p2, type = next.type(trig.order.id)));
		C2 = sum.pm(C2, mult.pm(pT1, dp.pm(p2, xn=xn)));
	}
	rlst = list(C1=C1, C2=C2, TrigLog=pT1, Trig0=pT0, Div=pT1);
	attr(rlst, "trig.order") = trig.order;
	attr(rlst, "D") = "TrigLog";
	return(rlst);
}
extract2.pm = function(p, s.warn=NULL) {
	check.f = function(p, sn)
		if(length(p) > 2) warning(paste0("Only the first 2 polynomials of ", sn, " processed!"));
	if(inherits(p, "pm")) {
		return(list(p, NULL));
	}
	if(is.list(pT)) {
		if(length(pT) >= 2) {
			p1 = p[[1]]; p2 = p[[2]];
			if( ! is.null(s.warn)) check.f(p, s.warn);
			return(list(p1, p2));
		} else return(list(pT[[1]], NULL));
	}
	# TODO: error;
}


################
### D( ODE ) ###
################

# D( dny ) => d(n+1)y;
dy.pm = function(p, yn="y", xn="x", starts.with=FALSE) {
	regEnd = if(starts.with) "" else "$";
	reg = paste0("^(?:d[0-9]*+|)", yn, regEnd);
	ynms = grepl(reg, names(p), perl=TRUE);
	ynms = names(p)[ynms];
	if(length(ynms) == 0) {
		# Note: NOT 0 ???
		if(is.null(xn)) return(p);
		return(dp.pm(p, xn=xn));
	}
	# Dx
	r = if(is.null(xn)) 0 else dp.pm(p, xn);
	# Dy
	dy.f = function(yn, dyn) {
		dp = p[p[,yn] != 0, , drop=FALSE];
		if(nrow(dp) == 0) return(0);
		dp$coeff = dp$coeff * dp[,yn];
		dp[,yn]  = dp[,yn] - 1;
		if( ! dyn %in% names(dp)) dp[,dyn] = 0;
		dp[,dyn] = dp[,dyn] + 1;
		return(dp);
	}
	y0 = if(starts.with) {
		ynms[startsWith(ynms, yn)];
	} else {
		ynms[ynms %in% yn];
	}
	if(length(y0) > 0) {
		for(y00 in y0) {
			dy = paste0("d", y00);
			r = sum.pm(r, dy.f(y00, dy));
		}
		ynms = ynms[ ! ynms %in% y0];
		if(length(ynms) == 0) return(r);
	}
	#
	dy = Dy.names(ynms);
	for(id in seq_along(ynms)) {
		r = sum.pm(r, dy.f(ynms[[id]], dy[[id]]));
	}
	#
	return(r);
}
Dy.names = function(nms) {
	# processes only "dy" & "dny";
	# dy.pm: "(?=.)" should NOT be necessary!
	npos = regexpr("^d[0-9]*+(?=.)", nms, perl=TRUE);
	len = attr(npos, "match.length");
	if(any(len < 0)) stop("Invalid dy name!");
	#
	n = rep(1, length(len));
	isDn = (len > 1);
	s = substr(nms[isDn], 2, len[isDn]); # pmax(2, len[isDn])
	n[isDn] = as.integer(s);
	n = n + 1;
	nms = substr(nms, len + 1, nchar(nms));
	nms = paste0("d", n, nms);
	return(nms)
}


#######################
#######################

# Split p1/p2 = p0 + D(p2)/p2 + ct/p2;
split.pm.fraction = function(p1, p2, by="x", drop=TRUE) {
	xn = by[1];
	idx  = match(xn, names(p1));
	idx2 = match(xn, names(p2));
	drop.f = function(p) if( ! drop) p else drop.pm(p);
	if(is.na(idx2) || (pow.max2 <- max(p2[, idx2])) == 0) {
		warning(paste0("p2 does NOT contain variable ", xn));
		return(list(P0=p1, D=0, Ct=0, Div.P0=p2));
	}
	pow.max = if(is.na(idx)) 0 else max(p1[, idx]);
	if(pow.max == 0) return(list(P0=0, D=0, Ct=p1, Div=p2));
	### Div
	pR = list(Div = p2);
	if(pow.max >= pow.max2) {
		tmp = div.pm(p1, p2, by=by); # may contain additional vars;
		if(nrow(tmp$Rem) == 0) {
			return(list(P0=tmp$Rez, D=0, Ct=0, Div=p2));
		}
		p1 = tmp$Rem; pR$P0 = drop.f(tmp$Rez);
	} else {
		pR$P0 = 0;
	}
	### D
	dp2 = dp.pm(p2, xn=xn);
	pow.Dmax = max(dp2[, xn]);
	if(pow.Dmax == 0) {
		pR$D = 0; pR$Ct = drop.f(p1); 
		return(pR);
	}
	tmp = div.pm(p1, dp2, by=by); # may contain additional vars;
	pR$D = drop.f(tmp$Rez);
	if(nrow(tmp$Rem) == 0) {
		pR$Ct = 0;
		return(pR);
	}
	pR$Ct = drop.f(tmp$Rem);
	return(pR);
}

### Solve Linear System

solve.LD.pm = function(pM, pR) {
	if(length(pM) > 4) stop("Not yet implemented!");
	pDiv = diff.pm(
		mult.pm(pM[[1]], pM[[4]]),
		mult.pm(pM[[2]], pM[[3]]));
	pC1 = diff.pm(
		mult.pm(pR[[1]], pM[[4]]),
		mult.pm(pM[[2]], pR[[2]]));
	pC2 = diff.pm(
		mult.pm(pM[[1]], pR[[2]]),
		mult.pm(pR[[1]], pM[[3]]));
	pR = list(C1=pC1, C2=pC2, Div=pDiv);
	return(pR);
}

####################
####################

### Integrate
I.pm = function(p, xn="x") {
	if( ! inherits(p, "data.frame")) stop("p must be a Polynomial!");
	idc = match("coeff", names(p));
	if(is.na(idc)) stop("p must be a Polynomial with coefficients!");
	# Variable:
	idx = match(xn, names(p));
	if(is.na(idx)) {
		warning(paste0("p does not contain the variable: ", xn));
		p[, xn] = 0;
		idx = match(xn, names(p));
	}
	# I
	p[, idx] = p[, idx] + 1;
	isLog = (p[, idx] == 0);
	p[ ! isLog, idc] = p[ ! isLog, idc] / p[ ! isLog, idx];
	if(any(isLog)) {
		warning("Logarithms generated!");
		nmLog = paste0("Log.", xn);
		p[, nmLog] = 0;
		p[isLog, nmLog] = 1;
	}
	return(p);
}
