
########################
###
### Leonard Mada
### [the one and only]
###
### P6 Polynomials
### Derived from P4
###
### draft v.0.1a


### Generate P6
### by Entangling roots of P4

### TODO:
# - compute elementary polynomials;


#####################

###############
### History ###

# draft v.0.1a:
# - moved from Polynomials.Derived.P6.R [v.0.3e [pre-z]]
#   to separate file;



####################

library(polynom)
library(pracma)

### helper functions
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
round0.p = function(p, tol=1E-7) {
	p = round0(as.vector(p), tol=tol)
	class(p) = "polynomial"
	return(p)
}
round0 = function(m, tol=1E-7) {
	m[abs(Re(m)) < tol & abs(Im(m)) < tol] = 0
	isZero = (Re(m) != 0) & (abs(Re(m)) < tol)
	if(sum(isZero) > 0) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	isZero = (Im(m) != 0) & (abs(Im(m)) < tol)
	if(sum(isZero) > 0) {
		m[isZero] = Re(m[isZero])
	}
	return(m)
}

###################

### P4 => P6

###
r4 = roots(c(1,0,0,1, 1))
r.g = expand.grid(r4, r4)
x = round0((r.g[,1] - r.g[,2])^2)
x = unique(x[ x != 0 ])
round0.p(poly.calc(x))
229 + 216*x - 112*x^2 + 26*x^3 + 8*x^4 + x^6

# TODO:
# x^6 + (-3*E1^2 + 8*E2)*x^5 + (...)*x^4 + ...
#
# b4 = (3*r^4 - 4*r[i]*r[j]^3 + 8*r[i]^2*r[j]^2 - 2*r[i]*r[j]*r[k]^2 + 12*E4)
#
# (3+4)*(E1^4 - 4*E1^2*E2 + 4*E1*E3 + 2*E2^2 - 4*E4)
# - 4*E1*(E1^3 - 3*E1*E2 + 3*E3) +
# + 8*r[i]^2*r[j]^2 + ...

# r^4 = E1^4 - 4*E1^2*E2 + 4*E1*E3 + 2*E2^2 - 4*E4
# r^3 = -(E1^3 - 3*E1*E2 + 3*E3)
# r^2 = E1^2 - 2*E2
#
# pow = 3
# sum(r4[1]*r4[-1]^pow, r4[2]*r4[-2]^pow, r4[3]*r4[-3]^pow, r4[4]*r4[-4]^pow)
#
# b3 = (-r^6 + 2*r[i]*r[j]^5 - 7*r[i]^2*r[j]^4 + 8*r[i]^3*r[j]^3 + 2*r[i]*r[j]^2*r[k]^3 + 16*r[ijk]*r[m]^3 +
#  -10*r[i]^2*r[j]^2*r[k]^2)*x^3 + ...

# Derivation:
# r^4 = E1*(E1^3 - 3*E1*E2 + 3*E3) - E2*(E1^2 - 2*E2) + E1*E3 - 4*E4
# (x-(r1-r2)^2)*(x-(r1-r3)^2)*(x-(r1-r4)^2)*(x-(r2-r3)^2)*(x-(r2-r4)^2)*(x-(r3-r4)^2)
# (x-(a-b)^2)*(x-(a-c)^2)*(x-(a-d)^2)*(x-(b-c)^2)*(x-(b-d)^2)*(x-(c-d)^2)


###
r4 = roots(c(1,2,0,0, 1))
r.g = expand.grid(r4, r4)
round0.p(poly.calc(r.g[,1] - r.g[,2]))
x = round0((r.g[,1] - r.g[,2])^2)
x = unique(x[ x != 0 ])
-176 + 288*x - 16*x^2 - 88*x^3 + 56*x^4 - 12*x^5 + x^6


###
r4 = roots(c(1,0,1,1, 2))
r.g = expand.grid(r4, r4)
round0.p(poly.calc(r.g[,1] - r.g[,2]))
x = round0((r.g[,1] - r.g[,2])^2)
x = unique(x[ x != 0 ])
1825 - 250*x - 335*x^2 + 86*x^3 + 38*x^4 + 8*x^5 + x^6


### P12 !
coeffs = c(1,2,0,0, 1)
alpha = 2
r4 = roots(coeffs)
r.g = expand.grid(r4, r4)
round0.p(poly.calc(alpha*r.g[,1] - r.g[,2])) / polynomial(rev(coeffs))
x = round0((alpha*r.g[,1] - r.g[,2]) )
x = unique(x[ x != 0 ])
# a P12 when alpha != 1


#################
### Special Cases


### m5 => P[2]*P[4]
m5 = unity(5, all=T)
m5 = m5[-1]

x.r = m5^3 - 2*m5^2 - m5
r.gr = expand.grid(x.r, x.r)
x = round0((r.gr[,1] - r.gr[,2])^2)
x = x[ x != 0]
15125 + 152500*x + 396700*x^2 + 62250*x^3 + 3740*x^4 + 100*x^5 + x^6


### ()^(1/4) => P[2]*P[4]
m4 = unity(4, all=T)
K = 2
k = K^(1/4) * m4

x.r = k^3 + k^2 + 3*k
r.gr = expand.grid(x.r, x.r)
x = round0((r.gr[,1] - r.gr[,2])^2)
x = x[ x != 0]
-658409500 + 11069440*x - 4889024*x^2 - 349696*x^3 + 16112*x^4 - 224*x^5 + x^6


