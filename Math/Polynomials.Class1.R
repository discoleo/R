
### Major Theories in Polynomials

### Polynomials of Class 1
### Roots with Simple Radicals
###
### Leonard Mada
### 2018 - 2020


##############
### Theory ###
##############

### TODO:
# - full theory!

# Base root:
# r = s[n-1]*k^(n-1) + s[n-2]*k^(n-2) + ... + s2*k^2 + s1*k + s0
# where k = K^(1/n) and s[j] and K are parameters;

# All roots:
# r[j] = sum(s[id]*k[j]^id)
# where k[j] = k * m^j, m^n = 1;
# and index id goes from 0 to n-1;


#################
### Exercises ###
#################

### E.1.) x^5 - 25/2 * x^2 - 125/4
# - let P(x) = x^5 - 25/2 * x^2 - 125/4;
# - let r0 = s*k^4 - s*k^3 + k^2 + k, where k = K^(1/5)
# E.1.a.) Show that for K=2 and s=1/2, r0 is a root of P(x) = 0;
# E.1.b.) Show that r1 = s*(k*m)^4 - s*(k*m)^3 + (k*m)^2 + k*m,
#         is also a root of P(x) = 0, where m^5 = 1, m = root of unity;
# E.1.c.) Find the remaining roots of this polynomial;

K = 2
s = 1/2
#
k = K^(1/5)
x = s*k^4 - s*k^3 + k^2 + k
err = x^5 - 25/2 * x^2 - 125/4
err # small precision error
# Students are invited to provide the mathematical proof!


##################

### some Examples:
m = complex(re=cos(2*pi/5), im=sin(2*pi/5))
m = m^(0:4)

###
# Parameter: can be modified
K = 5
# Roots
k = m * K^(1/5)
x = k^4 + k^3 + k
x^5 - 5*K* x^3 - 5*(K^2 + K) * x^2 - 5*K^3 * x - K^4 - K^3 - K - 5*(K^3 - K^2)

# for K=2
K = 2 # fixed
k = m * K^(1/5)
x = k^4 + k^3 + k
x^5 - 10*x^3 - 30*x^2 - 40*x - 46

# for K=3
K = 3 # fixed
k = m * K^(1/5)
x = k^4 + k^3 + k
x^5 - 15*x^3 - 60*x^2 - 135*x - 201

###
K = 3 # can be modified
# test
k = K^(1/5)
x = k^4 - 3*k^2 - k
x
x^5 + 5*K*x^3 + 15*K*(K + 3)*x^2 + 5*K*(28*K - 3)*x - K^4 + 15*K^3 + 198*K^2 + K


###
K = 2
s = 1/2
#
k = m*K^(1/5)
x = s*k^4 - s*k^3 + k^2 + k
#
x^5 - 5*K*(K*s^2 + K*s^3 - s + 1)*x^2 + 5*K*(K^2*s^4 + K*s^3 + 3*K*s^2 - K*s - 1)*x +
- K - K^2 - 10*K^2*s - 10*K^2*s^2 + 10*K^3*s^3 - 10*K^3*s^4 + K^3*s^5 - K^4*s^5
# for K = 2, s = 1/2
# x^5 - 25/2 * x^2 - 125/4





#################

### TODO:
# - clean;
# - expand;

library(polynom)

# alternative function to compute coefficients of polynomial
elemPoly = function(x, start=0, adjustSign=TRUE) {
	coeff = sapply(c(start, seq_along(x)), function(n){
		sum(apply(combn(x, n), 2, prod))
	})
	if(adjustSign) {
		if(start %% 2 ==0) {
			adj = c(1,-1)
		} else {
			adj = c(-1,1)
		}
		len = length(coeff)
		adj = rep(adj, len/2)
		if(len %% 2 == 1) {
			adj = c(adj, adj[1])
		}
		coeff[abs(coeff) < 1E-10 ] = 0
		coeff = coeff * adj
	}
	return(coeff)
}

poly.radical.create = function(root_coeff, K, n=length(root_coeff), doRound=TRUE) {
	k = ifelse(K < 0,
			ifelse(n %% 2 == 0, complex(re=0, im=(-K)^(1/n)), -(-K)^(1/n) ),
			K^(1/n))
	len = length(root_coeff) - 1
	n_1 = n - 1
	# root rotation
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	### Roots
	x = sapply(0:n_1, function(i) sum(root_coeff * (k*m^i)^(0:len)) )
	
	### Polynomial
	# does give a warning with complex roots
	# poly = poly.calc(x)
	poly = elemPoly(x) # * (-1)^((n:0)+ifelse(n %% 2 == 0, 0, 1) )
	# round: dangerous, but this is a minimal polynomial with integer coeffs,
	# if K and coeff are integers
	if(doRound) {
		poly = round(poly)
	}
	poly.coeff = as.numeric(round(poly))
	poly.str = paste(poly.coeff, "*x^", n:0, sep="")
	poly.str = poly.str[poly.coeff != 0]
	
	poly.list = list(
		r = x,
		poly.coeff = poly,
		poly = paste(poly.str, collapse="+"),
		n = n
	)
	return(poly.list)
}

### Examples

# E.1.) root = k^4 - 2*k^2 - k, where k = K^(1/5)
K = 2
p = poly.radical.create(c(0, -1, -2, 0, 1), K, n=5)
p
# 114 + 160*x + 80*x^2 + 10*x^3 + x^5
err = sapply(1:p$n, function(id) sum(p$poly * p$r[id]^(p$n:0)) )
err
# with package polynom
p.p = polynomial(rev(p$poly))
p.p
predict(p.p, p$r)


# E.2.) root = k^5 - k^3 - 3*k^2 + k, where k = K^(1/6)
K = 2
p = poly.radical.create(c(0, 1, -3, -1, 0, 1), K, n=6)
p
# 3726 + 2268*x + 594*x^2 + 108*x^3 - 18*x^4 + x^6
err = sapply(p$r, function(root) sum(p$poly * root^(p$n:0)) )
err
# with package polynom
p.p = polynomial(rev(p$poly))
p.p
predict(p.p, p$r)


# E.3.) root = k^6 - k^3 - 3*k^2 + k, where k = K^(1/7)
### NOTE: overflows easily with loss of precision!
K = 2
p = poly.radical.create(c(0, 1, -2, -1, 0, 0, 1), K, n=7)
p
# 818 + 672*x - 126*x^2 + 308*x^3 + 98*x^4 - 14*x^5 + x^7
err = sapply(p$r, function(root) sum(p$poly * root^(p$n:0)) )
err
# with package polynom
p.p = polynomial(rev(p$poly))
p.p
predict(p.p, p$r)



