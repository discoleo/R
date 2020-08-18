
########################
###
### Leonard Mada
### [the one and only]
###
### P6 Polynomials:
### Known Polynomials
###
### draft v.0.1a


### Introduction
# - a small selection of interesting polynomials;
# - exact solutions (sometimes partial) are available;
# - many other variants are available in the files:
#  -- Polynomials.Derived.P6.fromP4.R;
#  -- Polynomials.Derived.P6.R;


######################

### Series: K - 6*x + 4*x^3 + x^6
# - generalized solution *NOT* yet known;

### 2 - 6*x + 4*x^3 + x^6
# Polynomials.Derived.P6.fromP4.R
# s = 0; p = p3der_test2.gen(c(1,s,0,2-s), s2=-1)

### 6 - 6*x + 4*x^3 + x^6
# Polynomials.Derived.P6.R
# coeff = c(2,  3, -1,  1,  0); p = solve.p3cd(coeff)

### 1 - 6*x + 4*x^3 + x^6   # Partial
# (x + 1)*(-1 + 5*x + 5*x^2 + x^3 + x^4 + x^5)
# TODO: solution for P5;



###

### 9 - 8*x + 2*x^3 + x^6
# Polynomials.Derived.P6.R
# r = 1 + c(1i,-1i)*sqrt(2); x = sapply(r, function(r) roots(c(1, -2+2*r, -4, 3-2*r)))



### K + 3*n^5*x - 5*n^3*x^3 + 3*n*x^5 + x^6
# Polynomials.Derived.P6.R [4 & 5 parameter versions]
K = 2; n = 3; k = (K+n^6)^(1/3) * m3.all; x = sapply(k, function(k) roots(c(1, n, k-n^2)))
K + 3*n^5*x - 5*n^3*x^3 + 3*n*x^5 + x^6 # root is scaled by n;
#
K = 2; k = (1/K+1)^(1/3) * m3.all; x = 1 / sapply(k, function(k) roots(c(1, 1, k-1)))
K + 3*K*x - 5*K*x^3 + 3*K*x^5 + x^6

