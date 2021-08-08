########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Polynomial Generators


### Polynomial Generators
# - generate various types of polynomials;


###############
### History ###
###############


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
extend.spm = function(p, n=2, vb="be", vR="R", vS="S", sort=TRUE) {
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

