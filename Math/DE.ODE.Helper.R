########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Helper Functions
###
### draft v.0.1b-def


### History

### draft v.0.1b - v.0.1b-def:
# - code improvement: any() instead of sum();
# - minor change to line.tan();
### draft v.0.1a:
# - moved helper functions to this new file;
# - old file:
#   DE.ODE.Polynomial.R;



###################

#################
### Functions ###

### Plot
# - plot the tangent in a point;

### Mathematical / Numerical:
# - compute root order n;
# - round epsilon to 0;

####################

library(pracma)

### helper functions

### Tangent
line.tan = function(x, col="red", dx=5, p=p, dp=dp, ...) {
	slope = dp(x, ...)
	x.max = ifelse( (abs(x) >= 1), dx*x, dx*(sign(x) - sign(x)^2 + 1) );
	isInf = abs(slope) == Inf
	x.max[isInf] = x[isInf]
	p.x = p(x, ...)
	lines(c(x, x.max), c(p.x, p.x + (x.max-x)*slope), col=col)
	return(slope)
}

### Mathematical

rootn = function(r, n) {
	ifelse( (Im(r) == 0 & Re(r) >= 0), r^(1/n), - (-r)^(1/n) )
}
### round()
round0 = function(m, tol=1E-7) {
	m[abs(Re(m)) < tol & abs(Im(m)) < tol] = 0
	isNotNA =  ! is.na(m)
	isZero = (Re(m) != 0) & (abs(Re(m)) < tol)
	if(any(isZero[isNotNA])) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	isZero = (Im(m) != 0) & (abs(Im(m)) < tol)
	if(any(isZero[isNotNA])) {
		m[isZero] = Re(m[isZero])
	}
	return(m)
}
