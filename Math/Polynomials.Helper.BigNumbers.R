########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions


#######################

library(polynom)
library(pracma)
library(gmp) # BigNumbers

### helper Functions

# - hack to work with BigNumbers;
# - workaround for: aggregate();

aggregate0.pm = function(p) {
	# hack for bigz numbers
	coeff = p$coeff;
	p$coeff = seq(nrow(p));
	aggr.bignum = function(id) {
		x.df = data.frame(coeff=0);
		s = sum(coeff[id]);
		x.df$coeff = s;
	}
	p.r = aggregate(coeff~., p, aggr.bignum);
	return(p.r);
}

#############

### TODO:
# - cleanup duplicate code;

div.bigpm = function(p1, p2, by="x", debug=TRUE) {
	# very simple division
	xn = by;
	idx2 = match(xn, names(p2));
	if(is.na(idx2)) stop(paste0("P2 must contain the variable: ", xn));
	idx1 = match(xn, names(p1));
	if(is.na(idx1)) stop(paste0("P1 must contain the variable: ", xn));
	if( ! is.data.frame(p2)) p2 = as.data.frame(p2);
	#
	xpow2 = max(p2[,idx2]);
	pDx = p2[p2[,idx2] == xpow2, ];
	idc2 = match("coeff", names(p2));
	idc1 = match("coeff", names(p1));
	c2 = pDx[,idc2];
	pRez = as.data.frame(array(0, c(0,2)));
	names(pRez) = c(xn, "coeff");
	#
	if(nrow(pDx) == 1) {
		idn = match(names(pDx)[-idc2], names(p1));
		print(idn);
		if(any(is.na(idn))) stop(paste0("No matching variables: ", names(pDx)[is.na(idn)]));
		while(TRUE) {
			if(nrow(p1) == 0) break;
			xpow1 = max(p1[,xn]);
			if(xpow1 < xpow2) break;
			px1 = p1[p1[,xn] == xpow1, ];
			for(nc in seq_along(idn)) {
				px1[, idn[nc]] = px1[, idn[nc]] - pDx[, nc];
			}
			px1[, idc1] = px1[, idc1] / c2;
			pRez = sum.bigpm(pRez, px1);
			tmp = mult.bigpm(px1, p2); tmp$coeff = - tmp$coeff;
			p1 = sum.bigpm(p1, tmp);
		}
	}
	if(debug) {
		if(nrow(p1) > 0) print("Not divisible!")
		else print("Divisible!");
	}
	return(list(Rez=pRez, Rem=p1));
}
