
### S6 from P[6]: Symmetries

# Some Experiments

# IF (x1, ..., x6) is a solution:
# - then all permutations are also solutions: = 6! permutations;
# - let s12 = x1 + x2; s34 = x3 + x4; s56 = x5 + x6;
# - let S = s12 + s34 + s56;
# - let p2s = s12*s34 + s34*s56 + s56*s12;
# - let p3s = s12*s34*s56;
# - similar with: p12 = x1*x2; ...
# - then p2s, p3s, ... are invariant under 6*8 permutations;
# => deg(P(p2s)) = 6! / (6*8) = 15;

# As 15 > 6: not directly useful for solving P[6];


####################

### Helper Functions

# library(pracma)

source("Polynomials.Helper.R")

### P[15]
gen.p15 = function(x, collapse = FALSE) {
	xn = lapply(seq(2,6), function(id) {
		s1 = x[1] + x[id];
		x  = x[c(-1, -id)];
		s2 = x[1] + x[2];
		s3 = x[3] + x[4];
		p2s = s1*(s2+s3) + s2*s3;
		s2 = x[1] + x[3];
		s3 = x[2] + x[4];
		p2s = c(p2s, s1*(s2+s3) + s2*s3);
		s2 = x[1] + x[4];
		s3 = x[2] + x[3];
		p2s = c(p2s, s1*(s2+s3) + s2*s3);
		if(collapse) {
			p2s = sum(p2s);
		}
		return(p2s)
	});
	xn = unlist(xn);
	return(xn)
}
### P[10]
gen.p10 = function(x) {
	S = sum(x);
	xn = lapply(seq(10), function(id) {
		if(id <= 4) {
			s1 = x[1] + x[2] + x[id + 2];
		} else if(id <= 7) {
			s1 = x[1] + x[3] + x[id - 1];
		} else if(id != 10) {
			s1 = x[1] + x[4] + x[id - 3];
		} else {
			s1 = x[1] + x[5] + x[6];
		}
		s2 = S - s1;
		p2s = s1*s2;
		return(p2s)
	});
	xn = unlist(xn);
	return(xn)
}

###################

### Examples:

###
b  = c(1,0,0,0,0,-1,-1)
x0 = roots(b)

err = x0^6 - x0 - 1
round0(err)

### P[15]
x = gen.p15(x0)
poly.calc0(x, digits=8)

x^15 - 42*x^12 + 21*x^10 + 453*x^9 - 288*x^7 - 1232*x^6 - 353*x^5 + 96*x^4 +
	- 1728*x^3 + 792*x^2 - 1296*x - 32

### P[10]
x = gen.p10(x0)
poly.calc0(x, digits=8)

x^10 + 66*x^7 + 123*x^5 + 129*x^4 - 66*x^2 + 64*x + 1


### Other
x = gen.p15(x0, collapse = TRUE)
# Note: roots are complex conjugates;
poly.calc0(x, digits=11)
