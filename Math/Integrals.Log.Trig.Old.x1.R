

### Helper Constants

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;

###################

# Note:
# - Correct (robust) formulas are in file:
#   Integrals.Log.Trig.R;


### [previous]
n = 12; k = 5;
integrate(\(x) x * log(sin(x)), 0, pi * k/n)
n2 = floor((n-1)/2); idn = seq(n2)/n; even = 1 - (n %% 2);
sn = sin(2*pi*k*idn); cs = cos(2*pi*k*idn);
pracma::zeta(3) * (1 - 1/n^3 - even * (-1)^k * 7/n^3)/4 - (pi*k/n)^2 * log(2)/2 +
	+ sum(cs * (pracma::psi(2, idn) + pracma::psi(2, 1 - idn))) / (8*n^3) +
	- sum(sn * (pracma::psi(1, idn) - pracma::psi(1, 1 - idn))) * pi * k / (2*n^3);


# Note: NO direct relation to Clausen function;


################
################

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


################

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


################

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


##################
### Generalisation

# Note: Initial variant (NOT robust);

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


### Various Steps:

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


### Helper: Series
x  = pi/7; # Test
iN = seq(20000);
log(sin(x)) # ==
- sum( cos(2*iN*x) / iN ) - log(2);
#
log(cos(x)) # ==
- sum( (-1)^iN * cos(2*iN*x) / iN ) - log(2);


### Old Derivations

### on [0, pi/6]
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


### on [0, pi/5]
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

