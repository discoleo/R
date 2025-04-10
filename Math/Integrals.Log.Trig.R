

### Integrals: Log(Trig)

# - this file covers primarily
#   integrals of type: Log( Trig );


### Helper

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;
# Note:
# Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)
# Catalan = (psi(1, 1/4) - psi(1, 3/4)) / 4^2


###################
###################

### I( log(sin(x)) ) on various intervals

# Maths 505: A MONSTER INTEGRAL!!!
# int from 0 to infty arctan(x)/(x(x+1)(x^2+1))
# https://www.youtube.com/watch?v=u-FKjn_83l8
# Note: see 2nd & 3rd intermediate integrals;


### on [0, pi/2]
integrate(\(x) log(sin(x)), 0, pi/2)
integrate(\(x) log(cos(x)), 0, pi/2)
- pi/2*log(2)
#
integrate(\(x) log(1 + cos(x)), 0, pi/2)
integrate(\(x) log(1 + sin(x)), 0, pi/2)
- pi/2*log(2) + 2*Catalan


### on [0, pi/4]
integrate(\(x) log(sin(x)), 0, pi/4)
- pi/4*log(2) - Catalan/2

###
integrate(\(x) log(cos(x)), 0, pi/4)
- pi/4*log(2) + Catalan/2


### Varia:
integrate(function(x) log(cos(x) + sin(x)), 0, pi/2)
- pi/4*log(2) + Catalan
#
integrate(function(x) log(cos(x) - sin(x)), 0, pi/4)
- pi/8*log(2) - Catalan/2
#
integrate(function(x) log(cos(x) + sin(x)), 0, pi/4)
- pi/8*log(2) + Catalan/2


### I( log(sin(x)) ) on [0, pi/8]
C8 = integrate(\(x) log(x) / (x^2 + 1), 0, tan(pi/8))$value
# Note: see the sum(psi());
# C8 = - Catalan/4 - 2*sum(psi);
#
integrate(\(x) log(sin(x)), 0, pi/8)
- pi/8*log(2) - Catalan/8 + C8/2;
integrate(\(x) - log(cos(x)), 0, pi/8)$value + # see next;
	- pi*log(2)/4 - Catalan/4;

### I( log(cos(x)) )
integrate(\(x) log(cos(x)), 0, pi/8)
- pi*log(2) / 8 - Catalan/8 - C8/2;
dd = 512; # n = 8; dd = 8*n^2;
id = 1:3; sn = sin(2*pi*id / 8); sg = c(1,-1,1);
- pi*log(2) / 8 + sum(
	+ sg*sn * (pracma::psi(1, id/16) - pracma::psi(1, 1 - id/16)) +
	- sg*sn * (pracma::psi(1, 1/2 - id/16) - pracma::psi(1, 1/2 + id/16)) ) / dd;


# Derivation:
n = 8; # n = 10; # any integer >= 2;
integrate(\(x) log(cos(x)), 0, pi/n)
n2 = 2*n; isEven = (n %% 2 == 0);
ns = if(isEven) (n-1) / 2 else n/2; id = seq(ns);
sn = sin(2*pi*id / n); ni = 0:1000;
sg = - (-1)^id; sg2 = if(isEven) sg else -sg;
- pi*log(2)/n + 1/2 * sum(sapply(id,
	\(id0) sg[id0]*sn[id0] / (n2*ni + id[id0])^2 +
		-  sg[id0]*sn[id0] / (n2*ni + n2 - id[id0])^2 +
		- sg2[id0]*sn[id0] / (n2*ni + n - id[id0])^2 +
		+ sg2[id0]*sn[id0] / (n2*ni + n + id[id0])^2) )


### on [0, pi/3]

### I( log(cos(x)) )
integrate(\(x) log(cos(x)), 0, pi/3)
- pi*log(2)/3 + sqrt(3)/(4*36) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) +
	+ pracma::psi(1, 1/3) - pracma::psi(1, 2/3))

###
integrate(\(x) log(sin(x)), 0, pi/6)
- pi*log(2)/6 - sqrt(3)/(4*36) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) +
	+ pracma::psi(1, 1/3) - pracma::psi(1, 2/3))

### I( log(sin(x)) )
integrate(\(x) log(sin(x)), 0, pi/3)
- pi*log(2)/3 + sqrt(3)/(8*6^2) *
	( pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1, 7/12) +
	- 5 * pracma::psi(1, 1/6) + 5 * pracma::psi(1, 5/6) +
	- 3 * pracma::psi(1, 2/6) + 3 * pracma::psi(1, 4/6));

### I( log(tan(x)) )
integrate(\(x) log(tan(x)), 0, pi/3)
sin(pi/3)/(4*6^2) *
	( pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1, 7/12) 
	- 7 * pracma::psi(1, 1/6) + 7 * pracma::psi(1, 5/6) +
	- 5 * pracma::psi(1, 1/3) + 5 * pracma::psi(1, 2/3) );


# Derivation:
# for technique see from 12:40 in the Ref. Maths 505;
id = seq(0, 40000)
- pi*log(2)/3 + sqrt(3)/4 * sum(1/(6*id+1)^2, -1/(6*id+5)^2, 1/(6*id+2)^2, - 1/(6*id+4)^2)
- pi*log(2)/3 + sqrt(3)/(4*36) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) +
	+ pracma::psi(1, 2/6) - pracma::psi(1, 4/6))


### on [0, pi/5]
integrate(\(x) log(cos(x)), 0, pi/5)
sn = sin(2*pi/5 * c(1,2));
- pi*log(2)/5 +
	+ (sn[1] * (pracma::psi(1, 1/10) - pracma::psi(1, 9/10) +
		+ pracma::psi(1, 4/10) - pracma::psi(1, 6/10)) +
	- sn[2] * (pracma::psi(1, 2/10) - pracma::psi(1, 8/10) +
		+ pracma::psi(1, 3/10) - pracma::psi(1, 7/10))) / (2*100);

###
integrate(\(x) log(sin(x)), 0, pi/5)
integrate(\(x) - log(cos(x)), 0, pi * 3/10)$value - pi*log(2)/2;
sn = sin(6*pi/10 * c(1,2,3,4));
- pi*log(2)/2 - (- pi*log(2) * 3/10 +
	+ (sn[1] * (pracma::psi(1, 1/20) - pracma::psi(1, 19/20) +
		- pracma::psi(1, 9/20) + pracma::psi(1, 11/20)) +
	- sn[2] * (pracma::psi(1, 2/20) - pracma::psi(1, 18/20) +
		- pracma::psi(1, 8/20) + pracma::psi(1, 12/20)) +
	+ sn[3] * (pracma::psi(1, 3/20) - pracma::psi(1, 17/20) +
		- pracma::psi(1, 7/20) + pracma::psi(1, 13/20)) +
	- sn[4] * (pracma::psi(1, 4/20) - pracma::psi(1, 16/20) +
		- pracma::psi(1, 6/20) + pracma::psi(1, 14/20))
	) / (2*20^2));

# Note: also alternating signs;
- pi*log(2)/5 + 1/2 * sum(
	sn[1]/(10*id+1)^2, - sn[1]/(10*id+9)^2, - sn[2]/(10*id+2)^2, sn[2]/(10*id+8)^2,
	- sn[2]/(10*id+3)^2, sn[2]/(10*id+7)^2, sn[1]/(10*id+4)^2, - sn[1]/(10*id+6)^2)


### on [0, pi * 2/5]
integrate(\(x) log(cos(x)), 0, 2*pi/5)
sn = sin(2*pi/5 * c(2,4));
- pi*log(2) * 2/5 +
	+ (sn[1] * (pracma::psi(1, 1/10) - pracma::psi(1, 9/10) +
		+ pracma::psi(1, 4/10) - pracma::psi(1, 6/10)) +
	- sn[2] * (pracma::psi(1, 2/10) - pracma::psi(1, 8/10) +
		+ pracma::psi(1, 3/10) - pracma::psi(1, 7/10))) / (2*100);

###
integrate(\(x) log(sin(x)), 0, 2*pi/5)
integrate(\(x) - log(cos(x)), 0, pi/10)$value - pi*log(2)/2;
sn = sin(2*pi/10 * c(1,2,3,4));
- pi*log(2)/2 - (- pi*log(2) * 1/10 +
	+ (sn[1] * (pracma::psi(1, 1/20) - pracma::psi(1, 19/20) +
		- pracma::psi(1, 9/20) + pracma::psi(1, 11/20)) +
	- sn[2] * (pracma::psi(1, 2/20) - pracma::psi(1, 18/20) +
		- pracma::psi(1, 8/20) + pracma::psi(1, 12/20)) +
	+ sn[3] * (pracma::psi(1, 3/20) - pracma::psi(1, 17/20) +
		- pracma::psi(1, 7/20) + pracma::psi(1, 13/20)) +
	- sn[4] * (pracma::psi(1, 4/20) - pracma::psi(1, 16/20) +
		- pracma::psi(1, 6/20) + pracma::psi(1, 14/20))
	) / (2*20^2));

# Derivation:
- pi*log(2) * 2/5 + 1/2 * sum(
	sn[1]/(10*id+1)^2, - sn[1]/(10*id+9)^2, - sn[2]/(10*id+2)^2, sn[2]/(10*id+8)^2,
	- sn[2]/(10*id+3)^2, sn[2]/(10*id+7)^2, sn[1]/(10*id+4)^2, - sn[1]/(10*id+6)^2)


### on [0, pi/6]

### I( log(cos(x)) )
integrate(\(x) log(cos(x)), 0, pi/6)
- pi*log(2)/6 + sqrt(3)/(4*12^2) *
	( pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1,  7/12) +
	- pracma::psi(1, 1/6) + pracma::psi(1, 5/6) +
	+ pracma::psi(1, 2/6) - pracma::psi(1, 4/6));

### I( log(sin(x)) )
integrate(\(x) log(sin(x)), 0, pi/6)
- pi*log(2)/6 - sqrt(3)/(4*36) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) +
	+ pracma::psi(1, 1/3) - pracma::psi(1, 2/3));

### I( log(tan(x)) )
integrate(\(x) log(tan(x)), 0, pi/6)
- sin(pi/3)/(2*12^2) *
	( 3 * pracma::psi(1, 1/6) - 3 * pracma::psi(1, 5/6) +
	+ 5 * pracma::psi(1, 1/3) - 5 * pracma::psi(1, 2/3) +
	+ pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1,  7/12) );


#
id = seq(0, 40000)
sn = sin(2*pi/6 * c(1,2,3));
- pi*log(2)/6 + 1/2 * sum(
	sn[1]/(12*id+1)^2, - sn[1]/(12*id+11)^2, - sn[2]/(12*id+2)^2, sn[2]/(12*id+10)^2,
	sn[1]/(12*id+4)^2, - sn[1]/(12*id+8)^2, - sn[2]/(12*id+5)^2, + sn[2]/(12*id+7)^2)


### on [0, pi/12]

### I( log(cos(x)) )
integrate(\(x) log(cos(x)), 0, pi/12)
id = seq(5); sn = sin(2*id*pi/12); sg = - (-1)^id;
- pi*log(2)/12 + sum(sg * sn *
	(pracma::psi(1, id/24) - pracma::psi(1, 1 - id/24) +
	- pracma::psi(1, 1/2 - id/24) + pracma::psi(1, 1/2 + id/24)) ) / (2*24^2);

### I( log(sin(x)) )
integrate(\(x) log(sin(x)), 0, pi/12)
integrate(\(x) log(cos(x)), pi * 5/12, pi/2)
id = seq(5); sn = sin(2*5*id*pi/12); sg = - (-1)^id;
- pi*log(2)/12 - sum(sg * sn *
	(pracma::psi(1, id/24) - pracma::psi(1, 1 - id/24) +
	- pracma::psi(1, 1/2 - id/24) + pracma::psi(1, 1/2 + id/24)) ) / (2*24^2);

### I( log(tan(x)) )
integrate(\(x) log(tan(x)), 0, pi/12)
id = seq(5); sn = sin(2*id*pi/12) + sin(10*id*pi/12); sg = - (-1)^id;
- sum(sg * sn *
	(pracma::psi(1, id/24) - pracma::psi(1, 1 - id/24) +
	- pracma::psi(1, 1/2 - id/24) + pracma::psi(1, 1/2 + id/24)) ) / (2*24^2);


### on [0, 5/12 * pi]
integrate(\(x) log(cos(x)), 0, pi * 5/12)
id = seq(5); sn = sin(2*5*id*pi/12); sg = - (-1)^id;
- pi*log(2)*5/12 + sum(sg * sn * (
	+ pracma::psi(1, id/24) - pracma::psi(1, 1 - id/24) +
	- pracma::psi(1, 1/2 - id/24) + pracma::psi(1, 1/2 + id/24)) ) / (2*24^2);

###
integrate(\(x) log(sin(x)), 0, pi * 5/12)
integrate(\(x) - log(cos(x)), 0, pi/12)$value - pi/2*log(2);
id = seq(5); sn = sin(2*id*pi/12); sg = - (-1)^id;
- pi*log(2) * 5/12 - sum(sg * sn * (
	+ pracma::psi(1, id/24) - pracma::psi(1, 1 - id/24) +
	- pracma::psi(1, 1/2 - id/24) + pracma::psi(1, 1/2 + id/24)) ) / (2*24^2);


### on [0, pi/24]
integrate(\(x) log(cos(x)), 0, pi/24)
n = 24; id = seq((n-1)/2); sn = sin(2*id*pi/n); sg = - (-1)^id; dv = 2*n;
- pi*log(2)/n + sum(sg * sn *
	(pracma::psi(1, id/dv) - pracma::psi(1, 1 - id/dv) +
	- pracma::psi(1, 1/2 - id/dv) + pracma::psi(1, 1/2 + id/dv)) ) / (8*n^2);


### Ap / Cp Constants
# - useful for: I( log(gamma(x+k)) )
Ap = function(p, iter=10000) {
	n = iter; id = seq(n) + 1/p;
	sum(id * log(id)) - (n^2 + (2/p + 1)*n + 1/p^2 + 1/p + 1/6)*log(n + 1/p)/2 +
		+ n^2/4 + n/(2*p);
}
#
Ap(6) - Ap(-6)
log(6)/6 - sqrt(3)/(4*36*pi) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) +
	+ pracma::psi(1, 1/3) - pracma::psi(1, 2/3))


### Varia: Clausen Function
# see e.g.
# 1. Michael Penn: When "normal" trig functions aren't enough -- the Clausen function.
#    https://www.youtube.com/watch?v=5kN4oH8W1r8


############

### I( x^p * log(sin(x)) )

# Maths 505: ONE TOUGH INTEGRAL BOI: int 0 to pi/2 x ln(1+cos(x))
# https://www.youtube.com/watch?v=FUY4keujknc

### I( x * log(sin(x)) )
integrate(\(x) x * log(sin(x)), 0, pi/2)
- 1/8 * pi^2 * log(2) + 7/16 * pracma::zeta(3)

### I( x * log(cos(x)) )
integrate(\(x) x * log(cos(x)), 0, pi/2)
- 1/8 * pi^2 * log(2) - 7/16 * pracma::zeta(3)


### I( x^2 * log(sin(x)) )
integrate(\(x) x^2 * log(sin(x)), 0, pi/2)
- 1/24 * pi^3*log(2) + 3/16 * pi * pracma::zeta(3)

### I( x^2 * log(cos(x)) )
integrate(\(x) x^2 * log(cos(x)), 0, pi/2)
- 1/24 * pi^3*log(2) - 1/4 * pi * pracma::zeta(3);


### on [0, pi/4]

### I( x * log(sin(x)) )
integrate(\(x) x * log(sin(x)), 0, pi/4)
- 1/32 * pi^2 * log(2) - pi * Catalan / 8 + 35/128 * pracma::zeta(3)

### I( x * log(cos(x)) )
integrate(\(x) x * log(cos(x)), 0, pi/4)
- 1/32 * pi^2 * log(2) + pi * Catalan / 8 - 21/128 * pracma::zeta(3)

### I( x * log(tan(x)) )
integrate(\(x) x * log(tan(x)), 0, pi/4)
7/16 * pracma::zeta(3) - pi * Catalan / 4;


### I( x^2 * log(sin(x)) )
integrate(\(x) x^2 * log(sin(x)), 0, pi/4)
- pi^2 * (pi*log(2) + 6*Catalan) / (3*64) +
	+ 3/2^8 * pi*pracma::zeta(3) +
	+ (pracma::psi(3, 1/4) - pracma::psi(3, 3/4)) / 24 / 4^4;


### I( x^2 * log(cos(x)) )
integrate(\(x) x^2 * log(cos(x)), 0, pi/4)
integrate(\(x) 1/3 * x^3 * tan(x), 0, pi/4)$value +
	+ 1/3 * (pi/4)^3 * log(cos(pi/4));
- 1/(3*64) * pi^3*log(2) + pi^2*Catalan / 32 +
	+ 3/256 * pi * pracma::zeta(3) +
	- (pracma::psi(3, 1/4) - pracma::psi(3, 3/4)) / 24 / 4^4;
# Note:
# - see I( x^3 / cos(x) ) in file:
#   Integrals.Trig.Fractions.Poly.R;
# - see also below ( x^2 * log(tan(x)) ) & Wolframalpha;


### Helper

# on [pi/2, pi]
integrate(\(x) x * log(sin(x)), pi/2, pi)
integrate(\(x) (pi - x) * log(sin(x)), 0, pi/2)
-3/8 * pi^2*log(2) - 7/16 * pracma::zeta(3);

#
integrate(\(x) x^2 * log(sin(x)), pi/2, pi)
integrate(\(x) (pi - x)^2 * log(sin(x)), 0, pi/2)
- pi^2 * pi/2*log(2) +
	- 2*pi*(- 1/8*pi^2*log(2) + 7/16 * pracma::zeta(3)) +
	- 1/(3*8)*pi^3*log(2) + 3/16 * pi * pracma::zeta(3);
- 7/24 * pi^3*log(2) - 11/16 * pi * pracma::zeta(3);


# on [0, pi]
integrate(\(x) x * log(sin(x)), 0, pi)
- 1/2 * pi^2*log(2);

#
integrate(\(x) x^2 * log(sin(x)), 0, pi)
- 1/3*pi^3*log(2) - 1/2 * pi * pracma::zeta(3);


### on [0, pi/4]
integrate(\(x) x * log(cos(x)), 0, pi/4)
(pi*Catalan - 1/4 * pi^2 * log(2) - 21/16 * pracma::zeta(3)) / 8
#
integrate(\(x) x * log(sin(x)), 0, pi/4)
(35/16 * pracma::zeta(3) - 1/4 * pi^2 * log(2) - pi*Catalan) / 8

### I( x^2 * log(tan(x)) )
integrate(\(x) x^2 * log(tan(x)), 0, pi/4)
- ((pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128 +
	+ 3/2*pi^2 * Catalan) / 24;

# sum( log(sin(x)) on [0, pi/2], log(tan(x)) )
# =>

integrate(\(x) x^2 * log(sin(x)), 0, pi/4)
- (pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / (48*128) +
	- pi^2 * Catalan * 3/96 +
	+ (- 1/12 * pi^3*log(2) + 3/16 * pi * pracma::zeta(3)) / 16;


### on [0, pi/3]

### I( x * log(sin(x)) )
integrate(\(x) x * log(sin(x)), 0, pi/3)
cs = cos(2*pi/3); sn = sin(2*pi/3);
pracma::zeta(3)/4 - pracma::zeta(3) / (4*3^3) - (pi/3)^2 * log(2)/2 +
	+ sum(cs * (pracma::psi(2, 1/3) + pracma::psi(2, 1 - 1/3))) / (8*3^3) +
	- sum(sn * (pracma::psi(1, 1/3) - pracma::psi(1, 1 - 1/3))
		) * pi / (2*3^3);

### I( x * log(cos(x)) )
integrate(\(x) x * log(cos(x)), 0, pi/3)
integrate(\(x) (pi/2 - x) * log(sin(x)), pi/6, pi/2)
id = 1:2; ic = 1:3; sn = sin(2*pi*id/6); cs = cos(2*pi*ic/6);
- pi^2 * log(2) / 18 - 3/16 * pracma::zeta(3) +
	+ pi * sqrt(3)/(8*36) *
	( pracma::psi(1, 1/6) - pracma::psi(1, 5/6) +
	+ pracma::psi(1, 1/3) - pracma::psi(1, 2/3) ) +
	+ sum(cs * (pracma::psi(2, ic/6) - pracma::psi(2, 1/2 + ic/6))) / (8*6^3) +
	- sum(sn * (pracma::psi(1, id/6) - pracma::psi(1, 1 - id/6))
		) * pi / (2*6^3);

### I( x * log(tan(x)) )
integrate(\(x) x * log(tan(x)), 0, pi/3)
id = 1:2; ic = 1:3; sn = sin(2*pi*id/6); cs = cos(2*pi*ic/6);
pracma::zeta(3) * (7/16 - 1 / (4*3^3)) +
	+ sum(cos(2*pi/3) * (pracma::psi(2, 1/3) + pracma::psi(2, 1 - 1/3))) / (8*3^3) +
	- pi * sin(2*pi/3) * (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / (2*3^3) +
	- pi * sin(2*pi/3)/(4*36) *
	( pracma::psi(1, 1/6) - pracma::psi(1, 5/6) +
	+ pracma::psi(1, 1/3) - pracma::psi(1, 2/3) ) +
	- sum(cs * (pracma::psi(2, ic/6) - pracma::psi(2, 1/2 + ic/6))) / (8*6^3) +
	+ sum(sn * (pracma::psi(1, id/6) - pracma::psi(1, 1 - id/6))
		) * pi / (2*6^3);


### on [0, pi/6]

### I( x * log(sin(x)) )
integrate(\(x) x * log(sin(x)), 0, pi/6)
id = 1:2; ic = 1:3; sn = sin(2*pi*id/6); cs = cos(2*pi*ic/6);
pracma::zeta(3)/4 - (pi/6)^2 * log(2)/2 +
	+ sum(cs * (pracma::psi(2, ic/6) - pracma::psi(2, 1/2 + ic/6))) / (8*6^3) +
	- sum(sn * (pracma::psi(1, id/6) - pracma::psi(1, 1 - id/6))
		) * pi / (2*6^3);

### I( x * log(cos(x)) )
integrate(\(x) x * log(cos(x)), 0, pi/6)
id = 1:2; ic = 1:3; sn = sin(2*pi*id/6); cs = cos(2*pi*ic/6);
- pracma::zeta(3) * (3/16 + 1 / (16*3^3)) +
- (pi/6)^2 * log(2) / 2 +
	+ sum(cos(2*pi/3) * (pracma::psi(2, 1/3) + pracma::psi(2, 1 - 1/3)) / 32 +
		- sin(2*pi/3) * (pracma::psi(1, 1/3) - pracma::psi(1, 1 - 1/3)) * pi / 8
	) / 3^3 +
	- sum(cs * (pracma::psi(2, ic/6) - pracma::psi(2, 1/2 + ic/6))) / (8*6^3) +
	+ sum(sn * (pracma::psi(1, id/6) - pracma::psi(1, 1 - id/6))
		) * pi / (2*6^3);

### I( x * log(tan(x)) )
integrate(\(x) x * log(tan(x)), 0, pi/6)
ic = 1:3; id = 1:2; cs = cos(2*pi*ic/6); sn = sin(2*pi*id/6);
pracma::zeta(3) * (7/16 + 1 / (16*3^3)) +
	+ sum(cs * (pracma::psi(2, ic/6) - pracma::psi(2, 1/2 + ic/6))) / (4*6^3) +
	- sum(sn * (pracma::psi(1, id/6) - pracma::psi(1, 1 - id/6))) * pi / (6^3) +
	- sum(cos(2*pi/3) * (pracma::psi(2, 1/3) + pracma::psi(2, 2/3)) / 4 +
		- sin(2*pi/3) * (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) * pi
	) / (8*3^3);


### on [0, pi/5]

### I( x * log(sin(x)) )
integrate(\(x) x * log(sin(x)), 0, pi/5)
id = 1:2; idp = 2*id*pi/5; cs = cos(idp); sn = sin(idp);
pracma::zeta(3)/4 - (pi/5)^2 * log(2)/2 +
	+ pracma::psi(2, 1) / (8*5^3) +
	+ sum(cs * (pracma::psi(2, id/5) + pracma::psi(2, 1 - id/5))) / (8*5^3) +
	- sum(sn * (pracma::psi(1, id/5) - pracma::psi(1, 1 - id/5))
		) * pi / (2*5^3);

### I( x * log(cos(x)) )
integrate(\(x) x * log(cos(x)), 0, pi/5)
id = 1:2; idp = 2*id*pi/5;
cs = cos(idp); cs2 = cos(2*idp);
sn = sin(idp); sn2 = sin(2*idp);
- pracma::zeta(3) * (3/4 - 3/4 / 5^3) / 4 +
	- (pi/5)^2 * log(2)/2 +
	- sum((cs - cs2/4) * (pracma::psi(2, id/5) + pracma::psi(2, 1 - id/5))
		) / (8*5^3) +
	+ sum((sn - sn2/2) * (pracma::psi(1, id/5) - pracma::psi(1, 1 - id/5))
		) * pi / (2*5^3);

### I( x * log(tan(x)) )
integrate(\(x) x * log(tan(x)), 0, pi/5)
id = 1:2; idp = 2*id*pi/5;
cs = cos(idp); cs2 = cos(2*idp);
sn = sin(idp); sn2 = sin(2*idp);
pracma::zeta(3) * (1 + 3/4 - 3/4/5^3)/4 +
	+ pracma::psi(2, 1) / (8*5^3) +
	+ sum((cs - cs2/8) * (pracma::psi(2, id/5) + pracma::psi(2, 1 - id/5))
		) / (4*5^3) +
	- sum((sn - sn2/4) * (pracma::psi(1, id/5) - pracma::psi(1, 1 - id/5))
		) * pi / (5^3);


### Gen: on [0, pi/n]
n = 7 # ODD
integrate(\(x) x * log(tan(x)), 0, pi/n)
id = seq((n-1)/2); idp = 2*id*pi/n;
cs = cos(idp); cs2 = cos(2*idp);
sn = sin(idp); sn2 = sin(2*idp);
pracma::zeta(3) * (1 + 3/4 - 3/4/n^3)/4 +
	+ pracma::psi(2, 1) / (8*n^3) +
	+ sum((cs - cs2/8) * (pracma::psi(2, id/n) + pracma::psi(2, 1 - id/n))
		) / (4*n^3) +
	- sum((sn - sn2/4) * (pracma::psi(1, id/n) - pracma::psi(1, 1 - id/n))
		) * pi / (n^3);


### ODD Integer:
n = 7; k = 2; # k = ANY Integer;
integrate(\(x) x * log(sin(x)), 0, pi * k/n)
id = seq((n-1)/2); sn = sin(2*pi*k*id/n); cs = cos(2*pi*k*id/n);
pracma::zeta(3) * (1 - 1/n^3)/4 - (pi*k/n)^2 * log(2)/2 +
	+ sum(cs * (pracma::psi(2, id/n) + pracma::psi(2, 1 - id/n))) / (8*n^3) +
	- sum(sn * (pracma::psi(1, id/n) - pracma::psi(1, 1 - id/n))
		) * pi * k / (2*n^3);

### EVEN Integer:
n = 12;
integrate(\(x) x * log(sin(x)), 0, pi/n)
id = seq(n/2 - 1); sn = sin(2*pi*id/n); cs = cos(2*pi*id/n);
pracma::zeta(3) * (1/2 + 3/n^3)/2 - (pi/n)^2 * log(2)/2 +
	+ sum(cs * (pracma::psi(2, id/n) + pracma::psi(2, 1 - id/n))) / (8*n^3) +
	- sum(sn * (pracma::psi(1, id/n) - pracma::psi(1, 1 - id/n))
		) * pi / (2*n^3);

### Gen: on [0, pi * k/n]
n = 12; # EVEN integer;
k =  5;
integrate(\(x) x * log(sin(x)), 0, pi * k/n)
id = seq(n/2 - 1); sn = sin(2*pi*k*id/n); cs = cos(2*pi*k*id/n);
zf = if(k %% 2 == 0) -8 else 6;
pracma::zeta(3) * (1 + zf/n^3)/4 - (pi*k/n)^2 * log(2)/2 +
	+ sum(cs * (pracma::psi(2, id/n) + pracma::psi(2, 1 - id/n))) / (8*n^3) +
	- sum(sn * (pracma::psi(1, id/n) - pracma::psi(1, 1 - id/n))
		) * pi * k / (2*n^3);


### Derivation:

### from [0, pi/2]
integrate(\(x) x * log(tan(x)), pi/4, pi/2)
integrate(\(x) - x * log(cos(x)), 0, pi/2)$value +
integrate(\(x) pi/4 * log(cos(x+pi/4)), -pi/4, pi/4)$value +
integrate(\(x) pi/4 * log(tan(x)), pi/4, pi/2)$value;


# on [0, pi/4]
integrate(\(x) x^2 * (pi/2 - x) / cos(x), 0, pi/2)
(pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128 +
	+ 7/2*pi*pracma::zeta(3);
# by parts =>
integrate(\(x) 1/2 * (pi*x - 3*x^2) * log((1-sin(x))/(1+sin(x))), 0, pi/2)
integrate(\(x) -2*(12*x^2-4*pi*x+pi^2/4) * log(tan(x)), 0, pi/4)


### Helper
x  = pi/7; # Test
iN = seq(20000);
log(sin(x)) # ==
- sum( cos(2*iN*x) / iN ) - log(2);
#
log(cos(x)) # ==
- sum( (-1)^iN * cos(2*iN*x) / iN ) - log(2);


### Derivation: I( x * log(sin(x)) )

# on [0, pi/6]
iN = seq(10000);
integrate(\(x) x * log(sin(x)), 0, pi/6)
integrate(\(x) sapply(x, \(x) - sum( cos(2*iN*x) / iN )*x - log(2)*x), 0, pi/6)
integrate(\(x) sapply(x, \(x) sum( sin(2*iN*x) / iN^2 )) / 2, 0, pi/6)$value +
	- sum( sin(2*iN*pi/6) / iN^2 ) * pi/12 - (pi/6)^2 * log(2)/2;
sum( 1 / iN^3 ) / 4 - sum( cos(2*iN*pi/6) / iN^3 ) / 4 +
	- sum( sin(2*iN*pi/6) / iN^2 ) * pi/12 - (pi/6)^2 * log(2)/2;
id = 1:2; ic = 1:3;
sn = sin(2*pi*id/6); cs = cos(2*pi*ic/6);
pracma::zeta(3)/4 - (pi/6)^2 * log(2)/2 +
	+ sum(cs * (pracma::psi(2, ic/6) - pracma::psi(2, 1/2 + ic/6))) / (8*6^3) +
	- sum(sn * (pracma::psi(1, id/6) - pracma::psi(1, 1 - id/6))
		) * pi / (12*6^2);

# on [0, pi/3]
integrate(\(x) x * log(sin(x)), 0, pi/3)
integrate(\(x) sapply(x, \(x) - sum( cos(2*iN*x) / iN )*x - log(2)*x), 0, pi/3)
integrate(\(x) sapply(x, \(x) sum( sin(2*iN*x) / iN^2 )) / 2, 0, pi/3)$value +
	- sum( sin(2*iN*pi/3) / iN^2 ) * pi/6 - (pi/3)^2 * log(2)/2;
sum( 1 / iN^3 ) / 4 - sum( cos(2*iN*pi/3) / iN^3 ) / 4 +
	- sum( sin(2*iN*pi/3) / iN^2 ) * pi/6 - (pi/3)^2 * log(2)/2;
cs = cos(2*pi/3); sn = sin(2*pi/3);
pracma::zeta(3)/4 - (pi/3)^2 * log(2)/2 + pracma::psi(2, 1) / (8*3^3) +
	+ sum(cs * (pracma::psi(2, 1/3) + pracma::psi(2, 1 - 1/3))) / (8*3^3) +
	- sum(sn * (pracma::psi(1, 1/3) - pracma::psi(1, 1 - 1/3))
		) * pi / (2*3^3);

# Note:
((pracma::psi(2, 3/6) - pracma::psi(2, 1/2 + 3/6))) / (6^3) # ==
- pracma::zeta(3) / 18


# on [0, pi/5]
iN = seq(10000);
integrate(\(x) x * log(sin(x)), 0, pi/5)
sum( 1 / iN^3 ) / 4 - sum( cos(2*iN*pi/5) / iN^3 ) / 4 +
	- sum( sin(2*iN*pi/5) / iN^2 ) * pi/10 - (pi/5)^2 * log(2)/2;
id = c(1,2); cs = cos(2*id*pi/5); sn = sin(2*id*pi/5);
pracma::zeta(3)/4 - (pi/5)^2 * log(2)/2 + pracma::psi(2, 1) / (8*5^3) +
	+ sum(cs * (pracma::psi(2, id/5) + pracma::psi(2, 1 - id/5))) / (8*5^3) +
	- sum(sn * (pracma::psi(1, id/5) - pracma::psi(1, 1 - id/5))
		) * pi / (2*5^3);

#
iN = seq(10000); sg = rep(c(1,-1), length(iN)/2);
integrate(\(x) x * log(cos(x)), 0, pi/5)
integrate(\(x) sapply(x, \(x) sum( sg*cos(2*iN*x) / iN )*x - log(2)*x), 0, pi/5)
integrate(\(x) sapply(x, \(x) - sum( sg*sin(2*iN*x) / iN^2 )) / 2, 0, pi/5)$value +
	- sum( sin(4*iN*pi/5) / (2*iN)^2 ) * 2*pi/10 +
	+ sum( sin(2*iN*pi/5) / iN^2 ) * pi/10 - (pi/5)^2 * log(2)/2;
- pracma::zeta(3) * 3/4 / 4 + sum( cos(2*iN*pi/5) / iN^3 ) / 4 +
	- sum( cos(4*iN*pi/5) / iN^3 ) / 16 +
	- sum( sin(4*iN*pi/5) / iN^2 ) * pi/20 +
	+ sum( sin(2*iN*pi/5) / iN^2 ) * pi/10 - (pi/5)^2 * log(2)/2;
id = 1:2; idp = 2*id*pi/5;
cs = cos(idp); cs2 = cos(2*idp);
sn = sin(idp); sn2 = sin(2*idp);
- pracma::zeta(3) * 3/4 / 4 - (pi/5)^2 * log(2)/2 +
	+ (1 - 1/4) * pracma::zeta(3) / (4*5^3) +
	- sum((cs - cs2/4) * (pracma::psi(2, id/5) + pracma::psi(2, 1 - id/5))
		) / (8*5^3) +
	+ sum((sn - sn2/2) * (pracma::psi(1, id/5) - pracma::psi(1, 1 - id/5))
		) * pi / (2*5^3);


#####################

### Non-Trig Interval

### I( x * log(sin(x)) )
# Maths 505: More complex than I imagined
# https://www.youtube.com/watch?v=wUtFgYjis_0
# Note: convert sin(x) = sinh(1i*x) / 1i;

integrate(\(x) x * log(sin(x)), 0, 1)
integrate(\(x) x * Re(log(1i * exp(-1i*x) * (1 - exp(2i*x)) / 2)), 0, 1)
log(1i/2) / 2 - 1i/3 + integrate(\(x) x * Re(log(1 - exp(2i*x))), 0, 1)$value
log(1i/2) / 2 - 1i/3 + zeta(3)/4 + 1i/2 * Li2(exp(2i)) - Li3(exp(2i))/4

# TODO:
# - polylog2: implement polylog( COMPLEX );
# - alternative formula;


###################
###################

### I( log(|cos(x)|) / x^2 )
# Maths 505: 2 ridiculously awesome integrals!
# https://www.youtube.com/watch?v=SwizwPy-GmE

### I( log(|cos(x)|) / x^2 )
integrate(\(x) log(abs(cos(x))) / x^2, 0, Inf, subdivisions=8029, rel.tol=1e-5)
- pi/2

### I( log(|sin(x)|) / x^2 )
integrate(\(x) log(abs(sin(x))) / x^2 - log(x)/x^2, 0, Inf, subdivisions=4029, rel.tol=1e-5)
- pi/2


### Varia:
integrate(\(x) log(x)/x^2 * (1 - (x + 1)*exp(-x)), 0, Inf, subdivisions=8029, rel.tol=1e-5)
1 - Euler


# Varia:
integrate(\(x) - log(x) / (x^2 + x + 1), 0, 1)
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 9
# TODO: psi(1, x) - psi(1, 1-x) ???


#################

### I( log(tan(x)) )
# Michael Penn: Can you guess the trick for this integral?
# https://www.youtube.com/watch?v=8R0MiRYmjbk

# Intermediary:
integrate(function(x) log(tan(x)), 0, pi/4)
# ==
- Catalan


### I( log(cos x)^2 )
# Michael Penn: Integral of (ln(cos x))^2
# https://www.youtube.com/watch?v=ikyVHEHmgP8

integrate(\(x) log(cos(x))^2, 0, pi/2)
pi^3/24 + pi/2 * log(2)^2


#########################

### I( log(cos(x)) / tan(x) )
# Maths 505: A cool log trig integral
# https://www.youtube.com/watch?v=iPicVw4lOg0

# Note:
# Full Generalization at the end of this section;

### I( log(cos(x)) / tan(x) )
integrate(\(x) log(cos(x)) / tan(x), 0, pi/2)
integrate(\(x) 1/2 * log(1 - x^2) / x, 0, 1)
- pi^2/24

### I( log(cos(x)) / sin(x) )
integrate(\(x) log(cos(x)) / sin(x), 0, pi/2)
integrate(\(x) log(x) / (1 - x^2), 0, 1)
- pi^2/8

### I( log(cos(x)) / cos(x) )
integrate(\(x) log(cos(x)) / cos(x) - log(pi/2 - x)/(pi/2 - x), 0, pi/2)
integrate(\(x) pi/2*log(cos(pi/2*x)) / cos(pi/2*x) - (log(pi/2) + log(1-x))/(1 - x), 0, 1)
- pi^2/24 - log(pi)^2/ 2 + log(pi)*log(2);
- pi^2/24 + log(2)^2/2 - log(pi/2)^2/ 2;


# equivalent:
integrate(\(x) 1/2 * log(1 - x^2) / (1 - x^2) - 1/4*(log(1-x) + log(2)) / (1-x), 0, 1)
3*log(2)^2/8 - pi^2/24

### Variant: I( log(x^2+1) / sqrt(x^2+1) )
# tan(x) = y =>
# [is also directly computable: polynomial fraction]
integrate(\(x) log(x^2+1) / sqrt(x^2+1) - 2*log(x+1)/(x+1), 0, Inf)
- pi^2/12 + log(2)^2


# Lim: x -> 1
x = 1 - 1E-6;
log(1 - x^2) / (1 + x) - 1/2*(log(1-x) + log(2))
log1p(- x^2) / (1 + x) - 1/2*(log1p(-x) + log(2))
0;

# Derivation:
integrate(\(x) log(cos(x)) / cos(x) - log(pi/2 - x)/(pi/2 - x), 0, pi/2)
integrate(\(x) 1/4*(log(1-sin(x)) +
	+ log(2)) * cos(x) / (1-sin(x)) - log(pi/2 - x)/(pi/2 - x), 0, pi/2)$value +
	+ (3*log(2)^2/8 - pi^2/24);
- pi^2/24 + log(2)^2/2 - log(pi/2)^2/ 2;


### Generalisation

### Base:
p = -1/3; k = sqrt(5)
integrate(\(x) x^p / (x^2+1)^k, 0, Inf)
gamma((p+1)/2) * gamma(k - (p+1)/2) / gamma(k) / 2

### Base: Log
p = -1/3; k = sqrt(5)
integrate(\(x) x^p * log(x^2+1) / (x^2+1)^k, 0, Inf)
gamma((p+1)/2) * gamma(k - (p+1)/2) *
	(digamma(k) - digamma(k - (p+1)/2)) / gamma(k) / 2;


### I( sin(x)^p * cos(x)^q * log(cos(x)) )
p = -1/3; q = sqrt(5)
integrate(\(x) sin(x)^p * cos(x)^q * log(cos(x)), 0, pi/2)
gamma((p+1)/2) * gamma((q+1)/2) *
	(digamma((q+1)/2) - digamma((p+q+2)/2)) / gamma((p+q+2)/2) / 4

### I( sin(x)^p * cos(x)^q * log(sin(x)) )
p = -1/3; q = sqrt(5)
integrate(\(x) sin(x)^p * cos(x)^q * log(sin(x)), 0, pi/2)
gamma((p+1)/2) * gamma((q+1)/2) *
	(digamma((p+1)/2) - digamma((p+q+2)/2)) / gamma((p+q+2)/2) / 4


### I( sin(x)^p * cos(x)^q * log(cos(x))^2 )
p = -1/3; q = sqrt(5)
integrate(\(x) sin(x)^p * cos(x)^q * log(cos(x))^2, 0, pi/2)
gamma((p+1)/2) * gamma((q+1)/2) *
	(pracma::psi(1, (q+1)/2) - pracma::psi(1, (p+q+2)/2) +
	(digamma((q+1)/2) - digamma((p+q+2)/2))^2) / gamma((p+q+2)/2) / 8


### I( sin(x)^p * cos(x)^q * log(cos(x))^3 )
p = -1/3; q = sqrt(5)
integrate(\(x) sin(x)^p * cos(x)^q * log(cos(x))^3, 0, pi/2)
gamma((p+1)/2) * gamma((q+1)/2) * (
	pracma::psi(2, (q+1)/2)/2 - pracma::psi(2, (p+q+2)/2)/2 +
	+ (digamma((q+1)/2) - digamma((p+q+2)/2))/2 *
	(pracma::psi(1, (q+1)/2) - pracma::psi(1, (p+q+2)/2) +
		+ (digamma((q+1)/2) - digamma((p+q+2)/2))^2) +
	+ (digamma((q+1)/2) - digamma((p+q+2)/2)) *
		(pracma::psi(1, (q+1)/2) - pracma::psi(1, (p+q+2)/2))
	) / gamma((p+q+2)/2) / 8


### I( sin(x)^p * cos(x)^q * log(cos(x)) * log(sin(x)) )
p = -1/3; q = sqrt(5)
integrate(\(x) sin(x)^p * cos(x)^q * log(cos(x)) * log(sin(x)), 0, pi/2)
gamma((p+1)/2) * gamma((q+1)/2) *
	( - pracma::psi(1, (p+q+2)/2) +
	(digamma((p+1)/2) - digamma((p+q+2)/2)) *
	(digamma((q+1)/2) - digamma((p+q+2)/2)) ) / gamma((p+q+2)/2) / 8


### Special Cases:

### I( log(sin(x)) )
integrate(\(x) log(sin(x)), 0, pi/2)
integrate(\(x) log(cos(x)), 0, pi/2)
gamma(1/2)^2 * (digamma(1/2) + Euler) / 4
- pi*log(2)/2

### I( log(cos(x)) * log(sin(x)) )
integrate(\(x) log(cos(x)) * log(sin(x)), 0, pi/2)
- gamma(1/2)^2 * (pi^2/6 - 4*log(2)^2) / 8

### on [0, pi/4]
integrate(\(x) log(cos(x)) * log(sin(x)), 0, pi/4)
- gamma(1/2)^2 * (pi^2/6 - 4*log(2)^2) / 16


### I( cos(x)^2 * log(sin(x)) )
# Maths 505: This integral is INSANE beyond measure!
# https://www.youtube.com/watch?v=Em5R4ckyqk0

### I( cos(x)^2 * log(sin(x)) )
integrate(\(x) cos(x)^2 * log(sin(x)), 0, pi/2)
- pi*log(2)/4 - pi/8

### I( sin(x)^2 * log(sin(x)) )
integrate(\(x) sin(x)^2 * log(sin(x)), 0, pi/2)
gamma(1/2)^2 * (digamma(1/2) + Euler) / 4 + pi*log(2)/4 + pi/8;
- pi*log(2)/4 + pi/8;


### Log-Powers:

### I( log(cos(x))^2 )
integrate(\(x) log(cos(x))^2, 0, pi/2)
pi^3 / 24 + pi*log(2)^2 / 2


### I( log(cos(x))^3 )
integrate(\(x) log(cos(x))^3, 0, pi/2)
gamma(1/2)^2 * (
	pracma::psi(2, 1/2)/2 + pracma::zeta(3) +
	- log(2) * (3*pracma::psi(1, 1/2) - pi^2/2 + 4*log(2)^2)
	) / 8

### on [0, pi/4]
integrate(\(x) log(cos(x)) * log(sin(x)), 0, pi/4)
- gamma(1/2)^2 * (pi^2/6 - 4*log(2)^2) / 16
- pi^3 / 96 + pi*log(2)^2 / 4;


### Other: I( log(1 + sin(x)) * log(1 - sin(x)) )
# Maths 505: THE COOLEST LOG TRIG INTEGRAL ON YOUTUBE!
# https://www.youtube.com/watch?v=lS3yp-3QHRI

integrate(\(x) log(1 - sin(x)) * log(1 + sin(x)), 0, pi/2)
integrate(\(x) log(1 - sin(x)^2)^2 / 2, 0, pi/2)$value +
integrate(\(x) - log(1 - sin(x))^2 / 2 - log(1 + sin(x))^2 / 2, 0, pi/2)$value
- pi^3 / 12 + pi*log(2)^2 / 2;

### Helper
integrate(\(x) log(cos(x))^2 + log(sin(x))^2, 0, pi/4)
integrate(\(x) (log(cos(x))^2 + log(sin(x))^2) / 2, 0, pi/2)
integrate(\(x) log(cos(x))^2, 0, pi/2)
pi^3 / 24 + pi*log(2)^2 / 2

###
integrate(\(x) (log(1 - sin(x))^2 + log(1 + sin(x))^2) / 2, 0, pi/2)
integrate(\(x) log(2*cos(x)^2)^2 + log(2*sin(x)^2)^2, 0, pi/4)
pi^3 / 6 + pi*log(2)^2/2


# Weierstrass =>
integrate(\(x) 2*(2*log(1 + x)*log(1 - x) - log(1 + x^2)*log(1 - x^2)) / (x^2 + 1) +
	+ log(x^2 + 1)^2 / (x^2 + 1), 0, 1)
- pi^3/24 + pi/4 * log(2)^2

# Helper:
integrate(\(x) log(x^2 + 1)^2 / (x^2 + 1), 0, 1)
integrate(\(x) 4*log(cos(x))^2, 0, pi/4)
# TODO


#
integrate(\(x) log(sin(x))^2 - log(cos(x))^2, 0, pi/4)
integrate(\(x) log(sin(x)/(cos(x)+1)) * log(sin(x)/2) / 2, 0, pi/2)
integrate(\(x) log(tan(x)) * log(sin(2*x)), 0, pi/4)$value + Catalan*log(2);
integrate(\(x) - log(x) * log(x^2+1) / (x^2+1), 0, 1)$value + pi^3/16;
#
integrate(\(x) 1/2 * log(tan(x/2)) * log(sin(x)), 0, pi/2)
integrate(\(x) log(x) * log(2*x/(x^2+1)) / (x^2+1), 0, 1)
integrate(\(x) - log(x) * log(x^2+1) / (x^2+1), 0, 1)$value +
	+ log(2)*(pracma::psi(1,3/4) - pracma::psi(1,1/4)) / 16 +
	+ (pracma::psi(2,3/4) - pracma::psi(2,1/4)) / 64;
integrate(\(x) - log(x) * log(x^2+1) / (x^2+1), 0, 1)$value +
	- log(2)*Catalan + pi^3/16;

# TODO: ?


####################
####################

###
# Maths 505: Supreme integral solved using Feynman's trick (GONE WRONG!)
# https://www.youtube.com/watch?v=KpIplB9s9KM

### Helper:
# I( log(b^n + cos(x)^n) ) & Other

###
b = 7
integrate(\(x) log(b + cos(x)), 0, pi)
pi * acosh(b) - pi*log(2)

###
integrate(\(x) log(b^2 + cos(x)^2), 0, pi)
integrate(\(x) Re(log(b*1i + cos(x)) + log(- b*1i + cos(x))), 0, pi)
2*pi * Re(acosh(b * 1i)) - 2*pi*log(2)

### I( log(b^3 + cos(x)^3) )
b = sqrt(5)
integrate(\(x) log(b^2 + b*cos(x) + cos(x)^2), 0, pi)
2*pi * Re(acosh(b * (cos(pi/3) + 1i*sin(pi/3)))) - 2*pi*log(2)
#
integrate(\(x) log(b^3 + cos(x)^3), 0, pi)
2*pi * Re(acosh(b * (cos(pi/3) + 1i*sin(pi/3)))) +
	+ pi * acosh(b) - 3*pi*log(2)

### I( log(b^5 + cos(x)^5) )
b = sqrt(5)
integrate(\(x) log(b^5 + cos(x)^5), 0, pi)
2*pi * Re(acosh(b * (cos(pi/5) + 1i*sin(pi/5)))) +
2*pi * Re(acosh(b * (cos(3*pi/5) + 1i*sin(3*pi/5)))) +
	+ pi * acosh(b) - 5*pi*log(2)


### Compound:
integrate(\(x) log(1/cos(x)^2 + tan(x)^4), 0, pi/2)
pi * Im(acos(2 + 0i))
pi * log(sqrt(3) + 2)
pi/2 * acosh(7)

# Various forms of "same" Integral:
integrate(\(x) log(cos(x)^2 + sin(x)^4), 0, pi/2)
pi * log((sqrt(3) + 2) / 4)

integrate(\(x) log(1 + sin(x)^4/cos(x)^2), 0, pi/2)
pi * log((sqrt(3) + 2) / 2)

integrate(\(x) log(4 - sin(2*x)^2), 0, pi/2)
pi * log((sqrt(3) + 2) / 2)


### Polynomial
id = 2*pi * seq(5) / 11
r = 2*cos(id)
x = r^2 + 2*r - 1;
# poly.calc0(x)
x^5 - 2*x^4 - 16*x^3 - 9*x^2 + 26*x + 23
#
integrate(\(x) {
	x = cos(x);
	log(x^5 - 2*x^4 - 16*x^3 - 9*x^2 + 26*x + 23);
	}, 0, pi)
pi * sum(abs(Re(acosh(-x + 0i)))) +
	- 5*pi*log(2)


#########################

### I( sin(x)/x * log(1 + cos(x)) / cos(x) ) on [0, Inf]
# Maths 505: Feynman's technique is unreasonably OP!
# https://www.youtube.com/watch?v=xOtQ8Mh0cvg

# Note: upper = Inf
# but integral behaves very badly;
integrate(\(x) sin(x)/x * log(1 + cos(x)) / cos(x), 0, 200, subdivisions=4096, rel.tol=1E-5)
(pi/2 - log(2)) * pi/2


###
integrate(\(x) cos(x) / x - 1/x, 0, pi/2, rel.tol=1E-6)
- sum(pracma::expint(pi/2*1i*c(-1,1)))/2 - log(pi/2) - Euler

###
integrate(\(x) tan(x) / x - 2/pi / (pi/2 - x), 0, pi/2)
# TODO

# Note:
# lim( tan(x) / x - 2/pi / (pi/2 - x) ) = 4/pi^2
# for x -> pi/2;

#########################
#########################

### Derived Integrals

### I( log(x) / sqrt(2*b*x - x^2) )
# Maths 505: A nice integral from 1886
# https://www.youtube.com/watch?v=CbAcyU8gFPw

###
integrate(\(x) log(x) / sqrt(4*x - x^2), 0, 2)
- 2*Catalan

###
integrate(\(x) log(x) / sqrt(6*x - x^2), 0, 3)
pi/2 *  log(3/2) - 2*Catalan

###
b = sqrt(5)
integrate(\(x) log(x) / sqrt(2*b*x - x^2), 0, b)
pi/2 *  log(b/2) - 2*Catalan

# TODO: log(sin)-interval [0, pi/3];


#########################
#########################

### I( log(x^2 + log(b*cos(x))^2) )
# Maths 505: All my favourite advanced calculus tricks in one integral!
# https://www.youtube.com/watch?v=9cC2jKSclEY

integrate(\(x) log(x^2 + log(cos(x))^2), 0, pi/2)
pi * log(log(2))

###
b = 1/5
integrate(\(x) log(x^2 + log(b*cos(x))^2), 0, pi/2)
pi * log(log(2/b))


### Derivatives:

b = 1/5
integrate(\(x) log(b*cos(x)) / (x^2 + log(b*cos(x))^2), 0, pi/2)
- pi/2 / log(2/b)

# TODO:
# - explore more derivatives of integral;


#######################
#######################

### I( cos(2*x) / (log(tan(x)) * cos(x)^4 )
# Maths 505: A nasty looking trigonometric integral
# https://www.youtube.com/watch?v=wxSVkD9MQ-c
# Intermediate: I( (1 - z^2) / log(z) ) on [0, 1];

integrate(\(x) cos(2*x) / (log(tan(x)) * cos(x)^4), 0, pi/4)
- log(3)

