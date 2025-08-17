#####################
##
## Leonard Mada
## (the one and only)
##
## Beta Function / Distribution
## = Generalisations =
##
## v.0.1a


### Extensions & Generalisations


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;


####################

####################
### 2D Extension ###

# Note:
# - for various other variants with Radicals,
#   see file: Integrals.Double.Radicals.R;


### Gen Full: I( x^p * y^q / (1 - x^n*y^m)^k )
k = 1/sqrt(19);
m = sqrt(5); n = sqrt(7);
p = sqrt(3) - 1/5; q = sqrt(2) - 1/11;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p*y^q / (1 - x^m*y^n)^k, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta((p+1)/m, 1-k) - beta((q+1)/n, 1-k)) / ((q+1)*m-(p+1)*n);


### Reformulations:

### I( x^p * y^q / (1 - x^n*y^n)^(1/n) )
n = sqrt(19);
p = sqrt(3) - 1/5; q = sqrt(2) - 1/7;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p*y^q / (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta((p+1)/n, 1-1/n) - beta((q+1)/n, 1-1/n)) / ((q-p)*n);


### I( x^p * y^q * (1 - x^n*y^n)^(1/n) )
p = sqrt(2); q = -1/sqrt(3); n = 1/sqrt(7);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p*y^q * (1 - x^n*y^n)^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta((p+1)/n, 1+1/n) - beta((q+1)/n, 1+1/n)) / ((q-p)*n);


### I( x^p * y^q / (1 - x*y)^k )
k = 1/sqrt(19);
p = sqrt(3) - 1/5; q = sqrt(2) - 1/11;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^p*y^q / (1 - x*y)^k, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(beta(p+1, 1-k) - beta(q+1, 1-k)) / (q-p);

