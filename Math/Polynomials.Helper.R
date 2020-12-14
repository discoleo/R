########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions


#######################

library(polynom)
library(pracma)

### helper Functions

rootn = function(r, n) {
	if(n %% 2 == 0) return(r^(1/n)); # TODO: complex?
	ifelse( (Im(r) == 0 & Re(r) < 0), - (-r)^(1/n), r^(1/n) )
}
unity = function(n=3, all=TRUE) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	if(all) {
		m = m^(0:(n-1))
	}
	return(m)
}
mult.p = function(p1, p2) {
	p.m = outer(p1, p2)
    p = as.vector(tapply(p.m, row(p.m) + col(p.m), sum))
	return(p)
}
# round to 0
round0 = function(m, tol=1E-7) {
	m[abs(Re(m)) < tol & abs(Im(m)) < tol] = 0
	isZero = (Re(m) != 0) & (abs(Re(m)) < tol)
	if(any(isZero)) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	isZero = (Im(m) != 0) & (abs(Im(m)) < tol)
	if(any(isZero)) {
		m[isZero] = Re(m[isZero])
	}
	return(m)
}
round0.p = function(p, tol=1E-7) {
	p = round0(as.vector(p), tol=tol)
	class(p) = "polynomial"
	return(p)
}

solve.S = function(S, R, b=0) {
	# generic solver (based on existing S = x+y+z)
	b2 = if(length(b) > 1) b[2] else 0; # Ext A2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3;
	x = sapply(S, function(x) roots(c(1, -x, R[2] - b2*x, - R[3] + b3*x)))
	len = length(S)
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = R[3]/x - b3
	yz.s = S - x
	# TODO: robust (when necessary)
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

