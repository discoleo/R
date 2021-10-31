########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Generators
###
### draft v.0.1a


### ODE: Generators


####################

### Helper Functions


# include: Polynomials.Helper.R;
source("Polynomials.Helper.R")


### Test ODEs
# include: DE.ODE.Helper.R;
# source("DE.ODE.Helper.R")


#######################
#######################

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
	if(do.gcd) {
		xgcd = gcd.vpm(pD2R);
		if(xgcd > 1) pD2R$coeff = pD2R$coeff / xgcd;
	}
	pD2R = sort.dpm(pD2R, y="y");
	if(pD2R$coeff[1] < 0) pD2R$coeff = - pD2R$coeff;
	if(print) print(print.dpm(pD2R, do.sort=FALSE));
	return(pD2R);
}

# p1*sin(pT) + p2*cos(pT)
# where pT = pT0 + log(pT1)
genODE.TrigLog.pm = function(p1, p2, pT, f0=NULL, print=FALSE, pDiv=NULL, div.by=NULL) {
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
			d2f = dp.pm(df0, xn="x");
			if( ! isNZ.pm(d2f)) d2f = NULL;
		}
	}
	# convert Fractions from: D(log(...))
	pD2y = mult.pm(pR[[2]], pD$Div);
	pR[[2]] = pD2y;
	pR = solve.LD.pm(c(pC, pD[c("C1", "C2")]), pR);
	# lapply(pR, print.data.frame);
	# D2 =>
	pD$Div = NULL; # reset DIV; (could be useful in the future)
	pD2  = dp.trigLog.pm(pD);
	pD2R = mult.pm(pD2$C1, pR$C1); # pD2RC1
	pD2R = sum.pm(pD2R, mult.pm(pD2$C2, pR$C2)); # pD2RC2
	# d2y:
	pD2y = dy.pm(pD2y, yn="y", xn="x");
	if(! is.null(d2f)) {
		pD2y = diff.pm(pD2y, d2f);
	}
	pD2y = mult.pm(pD2y, mult.pm(pD2$Div, pR$Div));
	pD2R = diff.pm(pD2y, pD2R);
	if( ! is.null(pDiv)) pD2R = div.pm(pD2R, pDiv, by=div.by)$Rez;
	#
	pD2R = sort.dpm(pD2R, y="y");
	if(print) print(print.dpm(pD2R, do.sort=FALSE));
	return(pD2R);
}

