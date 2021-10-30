########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Derivation


### fast load:
# source("Polynomials.Helper.D.R")

### requires:
# - but is already loaded inside;
# source("Polynomials.Helper.R")

######################

###############
### History ###
###############


### draft v.0.1a:
# - moved Section on Derivation
#   from Polynomials.Helper.R;


#####################
#####################


# D( p(x) )
dp.pm = function(p, xn="x") {
	p = p[(p[,xn] != 0), , drop=FALSE];
	if(nrow(p) == 0) return(0);
	p$coeff = p$coeff * p[,xn];
	p[,xn] = p[,xn] - 1;
	return(p);
}
# D( p1 * exp(p2) )
dp.exp.pm = function(p, xn="x") {
	pr = mult.pm(dp.pm(p$Exp, xn), p$Poly);
	pr = sum.pm(pr, dp.pm(p$Poly, xn));
	p$Poly = pr;
	return(p);
}
# D( p1 / pdiv ])
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
			pT  = plst$TrigLog;
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
	pR = list(C1=pC1, C2=pC2, div=pDiv);
	return(pR);
}
