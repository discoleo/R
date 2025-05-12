

### Integrals: Log(Trig)

### Note:
# [refactor]
# - Old approaches moved to this file;
# - Approach based on Clausen function:
#   see file Integrals.Log.Trig.R;


### Helper

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;

#################

### on [0, pi/8]


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


################

### on [0, pi/3]

# Note:
# - see also section on [0, pi/6] (below);

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


################

### on [0, pi/5]

### I( log(cos(x)) )
integrate(\(x) log(cos(x)), 0, pi/5)
sn = sin(2*pi/5 * c(1,2));
- pi*log(2)/5 +
	+ (sn[1] * (pracma::psi(1, 1/10) - pracma::psi(1, 9/10) +
		+ pracma::psi(1, 4/10) - pracma::psi(1, 6/10)) +
	- sn[2] * (pracma::psi(1, 2/10) - pracma::psi(1, 8/10) +
		+ pracma::psi(1, 3/10) - pracma::psi(1, 7/10))) / (2*100);

### I( log(sin(x)) )
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


################

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


#################

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

