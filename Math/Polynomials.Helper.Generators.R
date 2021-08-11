########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Polynomial Generators
###
### draft v.0.1c-v2


### Polynomial Generators
# - generate various types of polynomials;


###############
### History ###
###############


### draft v.0.1c - v.0.1c-v2:
# - Generator for basic Class 2 polynomials;
### draft v.0.1b - v.0.1b-v2:
# - simple Generator for Class 1 polynomials;
### draft v.0.1a:
# - moved Generators from file:
#   Polynomials.Helper.R;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;


########################
########################

###############
### Class 1 ###
###############

### Simple Generator
# - roots = sum(b * k^seq(n-1)), where k^n = K;
# - accepts only numerical coefficients;
class1.simple.gen = function(b, n=length(b), kn="K", add.top=TRUE) {
	if(n < length(b)) warning("Invalid length of root coefficients");
	rename = function(p, name) {names(p)[1] = name; return(p);}
	nK = length(b);
	pK = data.frame(k=seq(0, nK-1), coeff=b);
	kn.lower = tolower(kn); # component of the root;
	pK = pK[pK$coeff != 0,];
	pK = rename(pK, kn.lower);
	pR = pK; lR = list();
	reduce.pow = function(p) {
		p[,kn] = p[,kn.lower] %/% n;
		p = p[, - match(kn.lower, names(p))];
		p = p[,c(2,1)];
		return(p)
	}
	for(pow in seq(n-1)) {
		pCoeff = pR[pR[, kn.lower] %% n == 0,];
		if(nrow(pCoeff) > 0) {
			pCoeff$coeff = -pCoeff$coeff * n / pow;
			pR = add.pm(pR, pCoeff);
			pCoeff = reduce.pow(pCoeff);
			lR[[n - pow]] = pCoeff;
		} else lR[[n - pow]] = 0;
		pR = mult.pm(pR, pK);
		pR = pR[pR$coeff !=0, ];
	}
	if(add.top) lR[[n]] = rename(data.frame(K=0, coeff=1), kn);
	pR = reduce.pow(pR);
	pR$coeff = - pR$coeff;
	lR$b0 = pR;
	attr(lR, "n") = n;
	return(lR);
}
toPoly.Class1S.pm = function(b, n=length(b), kn="K", xn="x") {
	p = class1.simple.gen(b, n=n, kn=kn, add.top=TRUE);
	lp = lapply(seq(length(p)), function(id) {
		if(is.data.frame(p[[id]])) cbind(p[[id]], x = if(id == (n+1)) 0 else id);
	});
	pR = do.call(rbind, lp);
	# pxtop = data.frame(x=n, coeff=1); names(pxtop)[1] = xn;
	# pR = sum.pm(pR, pxtop);
	return(pR);
}

### Class 2: Simple Type
# - from roots of unity of Order (n+1);
toPoly.Class2.pm = function(n, s.id=NULL, sn="s", xn="x", include.last=FALSE) {
	# include.last = if s[n] should be removed (as it is redundant);
	# s.id = offers greater control;
	len = if(include.last || ! is.null(s.id)) n else n-1;
	sn = paste0(sn, seq(0, len));
	ddf = as.data.frame(diag(len+1));
	ddf = rbind(ddf, 0);
	names(ddf) = sn;
	ddf$m = 0; ddf$m[seq(2, len+1)] = seq(len);
	ddf[,xn] = 0; ddf$coeff = -1;
	ddf[len+2, xn] = 1; ddf$coeff[len+2] = 1;
	if( ! is.null(s.id)) {
		if(is.logical(s.id)) {
			ddf = rbind(ddf[s.id,, drop=FALSE], ddf[len+2,, drop=FALSE]);
		} else if(is.numeric(s.id)) {
			ddf = ddf[c(s.id+1, len+2),];
		} else stop("Unsupported indexes!");
		ddf = reduce.var.pm(ddf);
	}
	#
	pR = ddf;
	for(pow in seq(2, n)) {
		p2 = ddf;
		p2$m = (pow * p2$m) %% (n+1);
		pR = mult.pm(pR, p2);
		pR$m = pR$m %% (n+1);
	}
	pR = aggregate0.pm(pR);
	# reduce sum of unity;
	pM = pR[pR$m == 1, ];
	idm = match("m", names(pM));
	for(nr in seq(nrow(pM))) {
		pm1 = pM[nr, -idm];
		pm1$coeff = - pm1$coeff;
		pm1 = as.data.frame(sapply(pm1, function(x) rep(x, n+1)));
		pm1$m = seq(0, n);
		pR = rbind(pR, pm1); # faster
		# pR = sum.pm(pR, pm1);
	}
	pR = aggregate0.pm(pR);
	pR = pR[pR$coeff != 0, ];
	if(all(pR$m == 0)) pR$m = NULL
	else print("Error: Roots of unity should have canceled out!")
	return(pR)
}



#######################
#######################

### Classic Polynomials
bR.gen = function(pb, pR=1) data.frame(b = pb, R = pR, coeff = 1)
bx.gen = function(pb, px=1) data.frame(b = pb, x = px, coeff = 1)
classic.BaseSimple.gen = function(n) data.frame(x=c(n, 1, 0), b=c(0,1,0), R=c(0,0,1), coeff=c(1,1,-1));
### S2 Simple Type: x^n + b*y = R
classic.S2Simple.gen = function(n=3) {
	p1 = data.frame(
		x = c(n,0), b = c(0,0), R = c(0,1),
		coeff = c(-1, 1)
	)
	p1 = pow.pm(p1, n);
	p1 = diff.pm(p1, bR.gen(n, pR=1));
	p1 = add.pm(p1, bx.gen(n+1, px=1));
	if(n %% 2 == 1) p1$coeff = - p1$coeff;
	p1 = sort.pm(p1, sort.coeff=c(4,2,3,1), xn="x")
	rownames(p1) = seq(nrow(p1))
	return(p1);
}

### Derived Polynomials
roots.derived = function(n, pow=seq(n-1), rn="r", sn="s", all.roots=TRUE) {
	slen = length(pow);
	S = diag(slen);
	s = lapply(seq(nrow(S)), function(nr) S[nr,]);
	s = data.frame(s); names(s) = paste0(sn, pow);
	p1 = data.frame(x=rep(0, slen), r1=pow, s, coeff=-1);
	names(p1)[2] = paste0(rn, 1);
	p1 = rbind(p1, c(1, rep(0, slen+1), 1))
	if(all.roots) {
		p.list = list(p1);
		for(id in seq(2, n)) {
			pT = p1;
			names(pT)[2] = paste0(rn, id);
			p.list[[id]] = pT;
		}
		p1 = p.list;
	}
	return(p1)
}

### Extensions
### Extensions to Eq S:
extend.spm = function(p, n=2, vb="be", vR="R", vS="S", sort=TRUE) {
	# R => (R - be[1]*S - be[2]*S^2);
	pS = data.frame(R=0, S=seq(n), coeff=-1);
	pS = rbind(pS, c(1,0,1));
	names(pS)[1:2] = c(vR, vS);
	b.all = paste0(vb, seq(n));
	for(id in seq(n)) {
		pB = data.frame(rep(0, n+1));
		pB[id,1] = id; names(pB) = b.all[id];
		pS = cbind(pS, pB);
	}
	# substitute in p
	p = replace.pm(p, pS, vR, pow=1);
	if(sort) p = sort.pm(p, c(4,3), xn=vS);
	return(p);
}



######################
######################

### Tests

### Class 1 Poly:
b = c(0,1,-1,0,2)
K = 2;
p = toPoly.Class1S.pm(b, 5)
r = sum((K^(1/5))^seq(4) * b[-1])

print.p(p, "x")
eval.pm(p, c(K, r));

