###################
##
## Trig: ATAN
## Polynomial Fractions of Atan
## Power = 4


####################

### Helper Constants
Catalan = 0.915965594177219015054603514;


####################
####################

### ATAN(...) / (x^4 + 1)

### on [0, Inf]

### I( atan(x) / (x^4 + 1) )
integrate(\(x) atan(x) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
(pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 64 +
	+ (digamma(5/8) - digamma(1/8)) * pi / 32;

### I( x^2 * atan(x) / (x^4 + 1) )
integrate(\(x) x^2 * atan(x) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
(pracma::psi(1, 3/8) - pracma::psi(1, 7/8)) / 64 +
	- (digamma(5/8) - digamma(1/8)) * pi / 32 + pi^2 / sin(3*pi/4) / 8;


### I( x * atan(k*x) / (x^4 + 1) )
k = 5^(1/3)
integrate(\(x) x * atan(k*x) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
cs = cos(pi/4); sn = sin(pi/4);
pi^2 / 8 - (atan((1/k + cs)/sn) - atan(cs/sn)) * pi/2;
pi^2 / 4 - pi/2 * atan((1/k + cs)/sn);

# Solution: Detailed formula
id = seq(1, 4, by=2) * pi/4; x = 1/k;
cs = cos(id); cs2 = cos(3*id); cs1 = cos(2*id);
sn = sin(id); sn2 = sin(3*id); sn1 = sin(2*id);
pi^2 / sin(2*pi/4) / 8 +
- (sum(cs2*log(x^2 + 2*cs*x + 1) +
	+ 2*sn2 * (atan((x + cs)/sn) - atan(cs/sn))) / sin(3*pi/4) +
sum(cs*log(x^2 + 2*cs*x + 1) +
	+ 2*sn * (atan((x + cs)/sn) - atan(cs/sn))) / sin(pi/4) +
sum(cs1*log(x^2 + 2*cs*x + 1) +
	+ 2*sn1 * (atan((x + cs)/sn) - atan(cs/sn))) * 2) * pi / 16;


### I( x^3 * atan(k*x) / (x^4 + 1) )
k = 3^(1/3); # k = - 1/5^(1/3);
integrate(\(x) x^3 * atan(k*x) / (x^4 + 1) - sign(k)*pi/2/(x+1), 0, Inf, rel.tol=1E-9)
(4*log(abs(k)) - log(k^4+1)) * sign(k) * pi/8 +
	+ log((k^2 - sqrt(2)*k + 1)/(k^2 + sqrt(2)*k + 1)) * pi/8;

# Detailed formula:
id = seq(1, 4, by=2) * pi/4; x = 1/k;
cs = cos(id); cs2 = cos(3*id);
sn = sin(id); sn2 = sin(3*id);
- log(x^4+1) * pi/8 +
(sum(cs2*log(x^2 + 2*cs*x + 1) +
	+ 2*sn2 * (atan((x + cs)/sn) - atan(cs/sn))) / sin(3*pi/4) +
- sum(cs*log(x^2 + 2*cs*x + 1) +
	+ 2*sn * (atan((x + cs)/sn) - atan(cs/sn))) / sin(pi/4) ) * pi/16;


###############

### ATAN-Pow: 2

### I( atan(x^2) / (x^4 + 1) )
integrate(\(x) atan(x^2) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
pi^2 / sin(pi/4) / 8 - pi * (digamma(1/2) - digamma(1/4)) / sin(pi/4) / 8
pi*(pi - 2*log(2)) * sqrt(2) / 16

### I( x^2 * atan(x^2) / (x^4 + 1) )
integrate(\(x) x^2 * atan(x^2) / (x^4 + 1), 0, Inf)
pi/8 / sin(pi/4) * (digamma(-1/4) - 2*digamma(-1/2) - Euler);


### Gen: I( x^p * atan(x^2) / (x^4 + 1) )
p = 1/sqrt(5); # p != 1;
integrate(\(x) x^p * atan(x^2) / (x^4 + 1), 0, Inf)
pi/8 / sin(pi*(p-1)/4) * (digamma((1-p)/4) - 2*digamma((1-p)/2) - Euler);

# Special Case:
integrate(\(x) atan(x^2) / (x * (x^4 + 1)), 0, Inf)
pi*log(2)/4;


### Gen: I( atan(k * x^2) / (x^4 + 1) )
k = 5^(1/3)
integrate(\(x) atan(k*x^2) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
x = 1/k^(1/2);
pi^2 / sin(pi/4) / 8 +
	+ (log(x^2 + 1) - 2*log(x + 1) - 2*atan(x)) * pi / sin(pi/4) / 8;

# Solution: Detailed formula
id = c(2,4) * pi/4; x = 1/k^(1/2);
cs = cos(id); cs1 = cs[1]; csp = cos(pi);
sn = sin(id); sn1 = sn[1]; snp = sin(pi);
pi^2 / sin(pi/4) / 8 +
- (sum(csp*log(x^2 - 2*cs1*x + 1) +
	- 2*snp * (atan((x - cs1)/sn1) + atan(cs1/sn1))) +
- (sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))))) * pi / sin(pi/4) / 8;


### Gen: I( x^2 * atan(k * x^2) / (x^4 + 1) )
k = 5^(1/3)
integrate(\(x) x^2 * atan(k * x^2) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
x = 1/k^(1/2);
pi^2 / sin(3*pi/4) / 8 +
	+ (2*log(x + 1) - log(x^2 + 1) - 2*atan(x)) * pi / sin(pi/4) / 8;

#  Solution: Detailed formula:
id = c(2,4) * pi/4; x = 1/k^(1/2);
cs = cos(id); cs1 = cs[1]; csp2 = cos(3*id); csp1 = cos(pi);
sn = sin(id); sn1 = sn[1]; snp2 = sin(3*id); snp1 = sin(pi);
pi^2 / sin(3*pi/4) / 8 +
- (sum(csp2*log(x^2 - 2*cs*x + 1) +
	- 2*snp2 * (atan((x - cs)/sn) + atan(cs/sn))) +
- (sum(csp1*log(x^2 - 2*cs1*x + 1) +
	- 2*snp1 * (atan((x - cs1)/sn1) + atan(cs1/sn1))))) * pi / sin(pi/4) / 8;


###############

### ATAN-Pow: 4

### I( atan(x^4) / (x^4 + 1) )
integrate(\(x) atan(x^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
pi^2 / sin(pi/4) / 8 - pi * (
	(digamma(13/16) - digamma(5/16)) / sin(pi/8) +
	(digamma(9/16) - digamma(1/16)) / sin(5*pi/8) +
	- (digamma(3/4) - digamma(1/4)) / sin(pi/4) * 2) / 32;


### I( x^2 * atan(x^4) / (x^4 + 1) )
integrate(\(x) x^2 * atan(x^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
pi^2 / sin(3*pi/4) / 8 - pi * (
	(digamma(15/16) - digamma(7/16)) / sin(3*pi/8) +
	(digamma(11/16) - digamma(3/16)) / sin(7*pi/8) +
	- (digamma(3/4) - digamma(1/4)) / sin(3*pi/4) * 2) / 32;


### Gen: I( atan(k * x^4) / (x^4 + 1) )
# see Integrals.Fractions.Unity.R;
k = 3^(1/4)
integrate(\(x) atan(k * x^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
# Solution:
pfr = seq(1, 8, by=2) * pi / 8; cs = cos(pfr); sn = sin(pfr);
cs4 = cos(5*pfr); sn4 = sin(5*pfr); x = 1 / k^(1/4);
cs3 = cos(4*pfr); sn3 = sin(4*pfr);
pi^2 / sin(pi/4) / 8 - pi*(
	+ sum(cs4*log(x^2 + 2*cs*x + 1) +
		+ 2*sn4 * (atan((x + cs)/sn) - atan(cs/sn))) / sin(pi/8) +
	+ sum(cs*log(x^2 + 2*cs*x + 1) +
		+ 2*sn * (atan((x + cs)/sn) - atan(cs/sn))) / sin(5*pi/8) +
	+ sum(cs3*log(x^2 + 2*cs*x + 1) +
		+ 2*sn3 * (atan((x + cs)/sn) - atan(cs/sn))) * 2 / sin(pi/4)) / 16;


### Gen: I( x^2 * atan(k * x^4) / (x^4 + 1) )
k = 3^(1/4)
integrate(\(x) x^2 * atan(k * x^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
# Solution:
pfr = seq(1, 8, by=2) * pi / 8; cs = cos(pfr); sn = sin(pfr);
cs6 = cos(7*pfr); sn6 = sin(7*pfr); x = 1 / k^(1/4);
cs3 = cos(4*pfr); sn3 = sin(4*pfr);
cs2 = cos(3*pfr); sn2 = sin(3*pfr);
pi^2 / sin(pi/4) / 8 - pi*(
	+ sum(cs6*log(x^2 + 2*cs*x + 1) +
		+ 2*sn6 * (atan((x + cs)/sn) - atan(cs/sn))) / sin(3*pi/8) +
	+ sum(cs2*log(x^2 + 2*cs*x + 1) +
		+ 2*sn2 * (atan((x + cs)/sn) - atan(cs/sn))) / sin(7*pi/8) +
	+ sum(cs3*log(x^2 + 2*cs*x + 1) +
		+ 2*sn3 * (atan((x + cs)/sn) - atan(cs/sn))) * 2 / sin(3*pi/4)) / 16;

# [alternative] Compact formula:
pi^2 / sin(pi/4) / 8 - pi*(
	sum((cs6/sin(3*pi/8) + cs2/sin(7*pi/8) + 2*cs3/sin(3*pi/4)) *
		log(x^2 + 2*cs*x + 1)) / 16 +
	sum((sn6/sin(3*pi/8) + sn2/sin(7*pi/8) + 2*sn3/sin(3*pi/4)) *
		(atan((x + cs)/sn) - atan(cs/sn))) / 8);


# Derivation:

# Base:
b = 5^(1/8)
integrate(\(x) 1/(x^2 + b^2) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) (1/(x^2 + b^2) - (x^2-b^2)/(x^4+1)) / (b^4 + 1), 0, Inf)
pi * (2/b + b^2/sin(pi/4) - 1/sin(3*pi/4)) / 4 / (b^4 + 1)

# from atan(x / b) / (x^4 + 1)
b = 5^(1/8)
integrate(\(x) x/(x^2 + b^2) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) (x/(x^2 + b^2) - (x^3 - b^2*x)/(x^4+1)) / (b^4 + 1), 0, Inf)
pi*b^2/sin(2*pi/4) / 4 / (b^4 + 1) - log(b) / (b^4 + 1)

# from x * atan(x / b) / (x^4 + 1)
b = 5^(1/8)
integrate(\(x) x^2/(x^2 + b^2) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) ((b^2*x^2 + 1)/(x^4+1) - b^2/(x^2 + b^2)) / (b^4 + 1), 0, Inf)
pi * (b^2/sin(3*pi/4) + 1/sin(pi/4) - 2*b) / 4 / (b^4 + 1)

# from x^2 * atan(x / b) / (x^4 + 1)
b = 5^(1/8)
integrate(\(x) x^3/(x^2 + b^2) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) ((b^2*x^3 + x)/(x^4+1) - b^2*x/(x^2 + b^2)) / (b^4 + 1), 0, Inf)
pi/sin(2*pi/4) / 4 / (b^4 + 1) + b^2 * log(b) / (b^4 + 1)

# from x^3 * atan(x / b) / (x^4 + 1)
b = 5^(1/8)
integrate(\(x) x^4/(x^2 + b^2) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) ((x^2 - b^2)/(x^4+1) + b^4/(x^2 + b^2)) / (b^4 + 1), 0, Inf)
pi * (1/sin(3*pi/4)- b^2/sin(pi/4) + 2*b^3) / 4 / (b^4 + 1)

# Helper:
integrate(\(b) log(b) / (b^4 + 1), 0, 1)
(pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64;


### Derivation: Higher Powers

# from atan(x^2 / b^2) / (x^4 + 1)
b = 5^(1/8) # does NOT include factor * 2*b
integrate(\(x) x^2/(x^4 + b^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) x^2 * (1/(x^4+1) - 1/(x^4 + b^4)) / (b^4 - 1), 0, Inf)
pi * (1/sin(3*pi/4) - 1/sin(pi/4) / b) / 4 / (b^4 - 1)

# from x^2 * atan(x^2 / b^2) / (x^4 + 1)
b = 5^(1/8) # does NOT include factor * 2*b
integrate(\(x) x^4/(x^4 + b^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) (b^4/(x^4 + b^4) - 1/(x^4+1)) / (b^4 - 1), 0, Inf)
pi * (b/sin(pi/4) - 1/sin(pi/4)) / 4 / (b^4 - 1)


# from atan(x^4 / b^4) / (x^4 + 1)
b = 5^(1/8) # does NOT include factor * 4*b^3
integrate(\(x) x^4/(x^8 + b^8) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) ((x^4+b^8)/(x^8 + b^8) - 1/(x^4+1)) / (b^8 + 1), 0, Inf)
pi * (1/sin(pi/8) * b + 1/sin(5*pi/8) / b^3 - 2/sin(pi/4)) / 8 / (b^8+1)

# from x^2 * atan(x^4 / b^4) / (x^4 + 1)
b = 5^(1/8) # does NOT include factor * 4*b^3
integrate(\(x) x^6/(x^8 + b^8) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) x^2 * ((x^4+b^8)/(x^8 + b^8) - 1/(x^4+1)) / (b^8 + 1), 0, Inf)
pi * (1/sin(3*pi/8) * b^3 + 1/sin(7*pi/8) / b - 2/sin(3*pi/4)) / 8 / (b^8+1)

# Note:
# (digamma(3/4) - digamma(1/4)) == pi;
# (allows slight simplification)


####################
####################

### on [0, 1]

### Pow = 4

### I( atan(x) / (x^4+1) )
integrate(\(x) atan(x) / (x^4+1), 0, 1)
integrate(\(x) -1/2 * log(x^2 + 1) / (x^4+1), 0, 1)$value +
(pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	+ (digamma(5/8) - digamma(1/8)) * pi / 16 +
	+ (digamma(5/8) - digamma(1/8)) * log(2) / 32 +
	- (digamma(3/4) - digamma(1/4)) *
		(digamma(7/8) - digamma(3/8)) / 64 +0;
# TODO


### I( x * atan(x) / (x^4+1) )
integrate(\(x) x * atan(x) / (x^4+1), 0, 1)
(digamma(7/8) - digamma(3/8)) *
	(digamma(5/8) - digamma(1/8)) / 64


### I( x^3 * atan(x) / (x^4+1) )
integrate(\(x) x^3 * atan(x) / (x^4+1), 0, 1)
(digamma(7/8) - digamma(3/8))^2 / 128 +
	- (digamma(5/8) - digamma(1/8))^2 / 128 + Catalan/2;


### Derivation:

# from atan(x/b) / (x^4 + 1)
b = sqrt(5)
integrate(\(x) x / (x^2+b^2) / (x^4+1), 0, 1)
integrate(\(x) x * (1/(x^2+b^2) - (x^2-b^2)/(x^4+1)) / (b^4+1), 0, 1)
integrate(\(x) (x/(x^2+b^2) - (x^3-b^2*x)/(x^4+1)) / (b^4+1), 0, 1)
(log(b^2 + 1)/2 - log(b) - log(2)/4 +
	+ b^2*(digamma(3/4) - digamma(1/4)) / 8) / (b^4+1)

