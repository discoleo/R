########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Generators
###
### draft v.0.1b


### ODE: Generators


####################

### Helper Functions


# include: Polynomials.Helper.R;
source("Polynomials.Helper.R")


### Test ODEs
# source("DE.ODE.Helper.R")


#######################
#######################

### Varia

eval.vpm = function(f, x) {
	sapply(x, function(x) eval.pm(f, x));
}


##################
### Generators ###
##################


### Trigonometric

# p1*sin(pT) + p2*cos(pT)
genODE.Trig.pm = function(p1, p2, pT, f0=NULL, pDiv=NULL, div.by=NULL,
		trig.order="sin", do.gcd=TRUE, print=FALSE) {
	if(is.null(p2)) p2 = data.frame(coeff=0);
	pC = list(p1, p2);
	pD = dp.trig.pm(pC, pT, trig.order=trig.order);
	# Linear System
	pR  = list(toPoly.pm("y"), toPoly.pm("dy"));
	d2f = NULL;
	if( ! is.null(f0)) {
		pR[[1]] = diff.pm(pR[[1]], f0);
		df0 = dp.pm(f0, xn="x");
		hasD = isNZ.pm(df0);
		if(hasD) {
			pR[[2]] = diff.pm(pR[[2]], df0);
			d2f = dp.pm(df0, xn="x");
			if( ! isNZ.pm(d2f)) d2f = NULL;
		}
	}
	# lapply(pR, print.pm);
	pR = solve.LD.pm(c(pC, pD[c("C1", "C2")]), pR);
	# D2 =>
	pD2 = dp.trig.pm(pD);
	pD2R = mult.pm(pD2$C1, pR$C1);
	pD2R = sum.pm(pD2R, mult.pm(pD2$C2, pR$C2));
	pD2y = pR$Div; pD2y$d2y = 1;
	if(! is.null(d2f)) {
		pD2y = diff.pm(pD2y, mult.pm(pR$Div, d2f));
	}
	pD2R = diff.pm(pD2y, pD2R);
	if( ! is.null(pDiv)) pD2R = div.pm(pD2R, pDiv, by=div.by)$Rez;
	#
	pD2R = format.dpm(pD2R, y="y", do.gcd=do.gcd);
	if(print) print(print.dpm(pD2R, do.sort=FALSE));
	return(pD2R);
}

# p1*sin(pT) + p2*cos(pT)
# where pT = pT0 + log(pT1)
genODE.TrigLog.pm = function(p1, p2, pT, f0=NULL, pDiv=NULL, div.by=NULL,
		do.gcd=TRUE, print=FALSE) {
	if(is.null(p2)) p2 = as.pm(0, x = "x", keep.zero = TRUE);
	if(is.numeric(p2)) p2 = as.pm(p2, x = "x", keep.zero = TRUE);
	pC = list(p1, p2);
	pD = dp.trigLog.pm(pC, pT);
	# Linear System
	pR  = list(toPoly.pm("y"), toPoly.pm("dy"));
	d2f = NULL;
	if( ! is.null(f0)) {
		pR[[1]] = diff.pm(pR[[1]], f0);
		df0 = dp.pm(f0, xn="x");
		hasD = isNZ.pm(df0);
		if(hasD) {
			pR[[2]] = diff.pm(pR[[2]], df0);
			# d2f = dp.pm(df0, xn="x");
			# if( ! isNZ.pm(d2f)) d2f = NULL;
		}
	}
	# convert Fractions from: D(log(...))
	pDy = mult.pm(pR[[2]], pD$Div);
	pR[[2]] = pDy;
	pR = solve.LD.pm(c(pC, pD[c("C1", "C2")]), pR);
	# lapply(pR, print.data.frame);
	# D2 =>
	pD$Div = NULL; # reset DIV; (could be useful in the future)
	pD2  = dp.trigLog.pm(pD);
	pD2R = mult.pm(pD2$C1, pR$C1); # pD2RC1
	pD2R = sum.pm(pD2R, mult.pm(pD2$C2, pR$C2)); # pD2RC2
	# d2y:
	pD2y = dy.pm(pDy, yn="y", xn="x");
	# if(! is.null(d2f)) {
		# pD2y = diff.pm(pD2y, d2f);
	# }
	pD2y = mult.pm(pD2y, mult.pm(pD2$Div, pR$Div));
	pD2R = diff.pm(pD2y, pD2R);
	if( ! is.null(pDiv)) pD2R = div.pm(pD2R, pDiv, by=div.by)$Rez;
	#
	powX = min(pD2R$x);
	if(powX != 0) pD2R$x = pD2R$x - powX;
	pD2R = as.pm(pD2R);
	pD2R = format.dpm(pD2R, y="y", do.gcd=do.gcd);
	if(print) print.dpm(pD2R, do.sort=FALSE);
	return(pD2R);
}

# p1*log(pL1) + p2*log(pL2)
genODE.Log.pm = function(p1, p2, pL1, pL2, f0=NULL, pDiv=NULL, div.by=NULL,
		do.gcd=TRUE, print=FALSE) {
	pC  = list(p1, p2);
	pC1 = list(C = p1, Log = pL1);
	pD1 = dp.log.pm(pC1);
	pC2 = list(C = p2, Log = pL2);
	pD2 = dp.log.pm(pC2);
	# convert Fractions from: D(log(...))
	pD  = expand.fr.pm(pD1, pD2);
	pDy = pD$Div; pDy$dy = 1;
	pDy = diff.pm(pDy, pD$B0);
	### Linear System
	py  = toPoly.pm("y");
	d2f = NULL;
	if( ! is.null(f0)) {
		py = diff.pm(py, f0);
		df0 = dp.pm(f0, xn="x");
		hasD = isNZ.pm(df0);
		if(hasD) {
			pDy = diff.pm(pDy, mult.pm(pD$Div, df0));
			d2f = dp.pm(df0, xn="x");
			if( ! isNZ.pm(d2f)) d2f = NULL;
		}
	}
	pR = list(py, pDy);
	pR = solve.LD.pm(c(pC, pD[c("C1", "C2")]), pR);
	# lapply(pR, print.data.frame);
	### D2 =>
	pD$Div = NULL; # reset DIV; (could be useful in the future)
	rename.f = function(l) { names(l) = c("C", "Log"); return(l); }
	pC1 = dp.log.pm(rename.f(pD[c("C1", "Log1")]));
	pC2 = dp.log.pm(rename.f(pD[c("C2", "Log2")]));
	pD2 = expand.fr.pm(pC1, pC2);
	# pD2 = (B0 + C1*Log1 + C2*Log2)/pD2$Div;
	pD2R = mult.pm(pD2$C1, pR$C1); # pD2RC1
	pD2R = sum.pm(pD2R, mult.pm(pD2$C2, pR$C2)); # pD2RC2
	pD2R = sum.pm(pD2R, mult.pm(pD2$B0, pR$Div));
	# d2y: d2f0 is already in pD2y;
	pD2y = dy.pm(pDy, yn="y", xn="x");
	pD2y = mult.pm(pD2y, mult.pm(pD2$Div, pR$Div));
	pD2R = diff.pm(pD2y, pD2R);
	if( ! is.null(pDiv)) pD2R = div.pm(pD2R, pDiv, by=div.by)$Rez;
	#
	pD2R = format.dpm(pD2R, y="y", do.gcd=do.gcd);
	if(print) print.dpm(pD2R, do.sort=FALSE);
	return(pD2R);
}

#################

### Helper Tools

format.dpm = function(p, y="y", do.gcd=TRUE) {
	if(do.gcd) {
		xgcd = gcd.vpm(p);
		if(xgcd > 1) p$coeff = p$coeff / xgcd;
	}
	p = sort.dpm(p, y=y);
	if(p$coeff[1] < 0) p$coeff = - p$coeff;
	return(p)
}

filter.names = function(l, exclude=NULL) {
	nms = names(l);
	if( ! is.null(exclude)) {
		isExcl = nms %in% exclude;
		nms = nms[ ! isExcl];
	}
	return(nms);
}

### Fractions
expand.fr.pm = function(p1, p2, add.names=TRUE) {
	# p[i] = (B0[i] + C[i]*FUN[i])/Div[i];
	pF1 = p1$Div; pF2 = p2$Div;
	# Free Polynomial:
	pS = mult.pm(pF2, p1$B0);
	pS = sum.pm(pS, mult.pm(pF1, p2$B0));
	# C * SomeFUN(...)
	pC1 = mult.pm(pF2, p1$C);
	pC2 = mult.pm(pF1, p2$C);
	pM  = mult.pm(p1$Div, p2$Div);
	pR  = list(B0=pS, C1=pC1, C2=pC2, Div=pM);
	# add _original_functions_;
	extract.other = function(p, id) {
		nms = filter.names(p, exclude=c("B0", "C", "Div"));
		lst = p[nms];
		names(lst) = paste0(names(lst), id);
		return(lst);
	}
	if(add.names) pR = c(pR, extract.other(p1, 1), extract.other(p2, 2));
	return(pR);
}

### Exponential Series
# exp(x) = sum( x^n / n! )
expand.Exp = function(n, asDiv=TRUE, xn="x") {
	div = factorial(n);
	p = data.frame(x = seq(n, 0));
	names(p) = xn;
	coeff = cumprod(c(1, rev(seq(n))));
	p$coeff = coeff;
	if( ! asDiv) { p$coeff = p$coeff / div; div = 1; }
	class(p) = c("pm", class(p));
	return(list(P=p, Div=div));
}

