

####################

### Helper Constants

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;

### Notes:
# - Function polylog2 is in file:
#   Integrals.Polylog.Helper.R;
# - Li2 is also available in package Rmpfr;


###################
###################

### Sums: 1 / (m^p * n^q * (m+n)^k )

# Sum( 1 / (m * n * (m+n)) )
id = 1:5000;
sum(sapply(id, \(id1) { sum(1 / (id*id1*(id+id1 + 0.0))); } ))
2*pracma::zeta(3);

# Sum( 1 / (m^2 * n * (m+n)) )
ids = 1:20000;
id1 = 1:5000; id2 = id1^2;
sum(sapply(ids, \(ids) { sum(1 / (id2*ids*(id1+ids))); } ))
integrate(\(x) sapply(x, \(x) - polylog2(x) * log(1-x) / x), 0, 1, rel.tol=1E-12)
pi^4 / 72;

# Sum( 1 / (m^2 * n^2 * (m+n)) )
id  = 1:2000;
id2 = id^2;
sum(sapply(id, \(id1) { 1 / (id1^2*id2*(id1+id)); } ))
integrate(\(x) sapply(x, \(x) polylog2(x)^2 / x), 0, 1, rel.tol=1E-12)
2*pracma::zeta(2) * pracma::zeta(3) - 3*pracma::zeta(5);

# Sum( 1 / (m^3 * n^3 * (m+n)) )
integrate(\(x) sapply(x, \(x) Re(polylog2(x, 3))^2 / x), 0, 1, rel.tol=1E-12)
4*pracma::zeta(7) - 2*pracma::zeta(2) * pracma::zeta(5);


# Sum( 1 / (m^2 * n * (m+n)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sapply(x, \(x) - polylog2(x*y) * log(1-x*y) / (x*y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pracma::zeta(2) * pracma::zeta(3) - 3/2 * pracma::zeta(5);

# Sum( 1 / (m^2 * n^2 * (m+n)^2) )
id  = 1:2000;
id2 = id^2;
sum(sapply(id, \(id1) { 1 / (id1^2*id2*(id1+id)^2); } ))
pracma::zeta(6) / 3;

# Integral-form:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sapply(x, \(x) polylog2(x*y)^2 / (x*y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pracma::zeta(6) / 3;


# Sum( 1 / (m^3 * n^2 * (m+n)) )
id  = 1:2000;
id2 = id^2;
sum(sapply(id, \(id1) { 1 / (id1^3*id2*(id1+id)); } ))
integrate(\(x) sapply(x, \(x) polylog2(x) * Re(polylog2(x, 3)) / x), 0, 1, rel.tol=1E-12)
pracma::zeta(3)^2 / 2;

# Sum( 1 / (m^3 * n^2 * (m+n)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sapply(x, \(x) Re(polylog2(x*y, 3)) * polylog2(x*y, 2) / (x*y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
2*pracma::zeta(7) - pracma::zeta(2) * pracma::zeta(5);


########################
########################

### Sums: 1 / (m^p * n^q * (m^2 + n^2) )

# Sum( 1 / (m^2 * n^2 * (m^2 + n^2)) )
id1 = (1:4000)^2;
id2 = id1;
zlp = sum(sapply(id1, \(id1) { 1 / (id1*id2*(id1+id2)); } ))
0.6543415088859439; # higher precision;
# TODO
