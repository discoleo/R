########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Asymmetric: Helper Functions
###
### draft v.0.1a


### specific Helper Functions
# - for Poly.System.Asym.Solv.S3.R;


################################
################################

### Order 4:
### A * (x^4+y^4+z^4) + B %*% c(x*y, x*z, y*z) = R


coeff.S4E2.S3P4 = function(r.m) {
	r  = r.m[,1]; r1 = r[1]; r2 = r[2]; r3 = r[3];
	ra = r.m[,2]; ra1 = ra[1]; ra2 = ra[2]; ra3 = ra[3];
	# helper function: simplistic implementation;
	coeff = c(
		- ((ra1*ra2)^4+(ra1*ra3)^4+(ra2*ra3)^4), # S^8
		(ra1*ra2*ra3)^2 +
			- 4*((ra1*ra2)^3*(ra1*r2+ra2*r1) + (ra1*ra3)^3*(ra1*r3+ra3*r1) +
			+ (ra2*ra3)^3*(ra2*r3+ra3*r2)), # S^7
	+ 2*(ra1*ra2*ra3)*(ra1*ra2*r3+ra1*ra3*r2+ra2*ra3*r1) +
		- 4*((ra1*ra2)^3*r1*r2 + (ra1*ra3)^3*r1*r3 + (ra2*ra3)^3*r2*r3) +
		- 6*((ra1*ra2)^2*(ra1*r2+ra2*r1)^2 + (ra1*ra3)^2*(ra1*r3+ra3*r1)^2 +
		+ (ra2*ra3)^2*(ra2*r3+ra3*r2)^2), # S^6
	+ (ra1*ra2*r3+ra1*ra3*r2+ra2*ra3*r1)^2 +
		+ 2*(ra1*ra2*ra3)*(ra1*r2*r3+ra2*r1*r3+ra3*r1*r2) +
		- 12*((ra1*ra2)^2*(ra1*r2+ra2*r1)*r1*r2 + (ra1*ra3)^2*(ra1*r3+ra3*r1)*r1*r3 +
			+ (ra2*ra3)^2*(ra2*r3+ra3*r2)*r2*r3) +
		- 4*((ra1*ra2)*(ra1*r2+ra2*r1)^3 + (ra1*ra3)*(ra1*r3+ra3*r1)^3 +
			+ (ra2*ra3)*(ra2*r3+ra3*r2)^3), # S^5
	- 6*((ra1*ra2)^2*(r1*r2)^2 + (ra1*ra3)^2*(r1*r3)^2 + (ra2*ra3)^2*(r2*r3)^2) +
		- 12*((ra1*ra2)*(ra1*r2+ra2*r1)^2*r1*r2 + (ra1*ra3)*(ra1*r3+ra3*r1)^2*r1*r3 +
			+ (ra2*ra3)*(ra2*r3+ra3*r2)^2*r2*r3) +
		- ((ra1*r2+ra2*r1)^4 + (ra1*r3+ra3*r1)^4 + (ra2*r3+ra3*r2)^4) +
		+ 2*(ra1*ra2*ra3)*(r1*r2*r3)
		+ 2*(ra1*ra2*r3+ra1*ra3*r2+ra2*ra3*r1)*(ra1*r2*r3+ra2*r1*r3+ra3*r1*r2), # S^4
	(ra1*r2*r3+ra2*r1*r3+ra3*r1*r2)^2 +
		+ 2*(ra1*ra2*r3+ra1*ra3*r2+ra2*ra3*r1)*(r1*r2*r3) +
		- 12*((ra1*ra2)*(ra1*r2+ra2*r1)*(r1*r2)^2 + (ra1*ra3)*(ra1*r3+ra3*r1)*(r1*r3)^2 +
			+ (ra2*ra3)*(ra2*r3+ra3*r2)*(r2*r3)^2) +
		- 4*((ra1*r2+ra2*r1)^3*r1*r2 + (ra1*r3+ra3*r1)^3*r1*r3 + (ra2*r3+ra3*r2)^3*r2*r3), # S^3
	2*(ra1*r2*r3+ra2*r1*r3+ra3*r1*r2)*(r1*r2*r3) +
		- 4*((ra1*ra2)*(r1*r2)^3 + (ra1*ra3)*(r1*r3)^3 + (ra2*ra3)*(r2*r3)^3) +
		- 6*((ra1*r2+ra2*r1)^2*(r1*r2)^2 + (ra1*r3+ra3*r1)^2*(r1*r3)^2 +
		+ (ra2*r3+ra3*r2)^2*(r2*r3)^2), # S^2
	(r1*r2*r3)^2 +
		- 4*((ra1*r2+ra2*r1)*(r1*r2)^3 + (ra1*r3+ra3*r1)*(r1*r3)^3 +
		+ (ra2*r3+ra3*r2)*(r2*r3)^3), # S^1
	- ((r1*r2)^4 + (r1*r3)^4 + (r2*r3)^4)
	)
	return(coeff);
}

### Test
R = c(1,2,3)
B = matrix(c(1,2,-1, 3,3,1, -1,2,-2), ncol=3, byrow=TRUE)
a = c(1,1,1)

x = -0.3481807989 - 1.1887367859i;
y =  0.6324344812 - 1.3543486953i;
z =  0.1302306191 + 1.1651128145i;
S4 = x^4 + y^4 + z^4;
m.coeff = solve(B, cbind(R, -1 * a));
r  = m.coeff[,1]; r1 = r[1]; r2 = r[2]; r3 = r[3];
ra = m.coeff[,2]; ra1 = ra[1]; ra2 = ra[2]; ra3 = ra[3];
