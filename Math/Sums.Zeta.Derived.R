


### Helper Functions

constEuler = 0.57721566490153286060651209008240243079;
Euler = constEuler;


###
zeta = function(n) {
	(1 + 1/(2^(n-1) - 1)) * pracma::eta(n);
}
unity.neg = function(n) {
	id = seq(1, 2*n, by=2);
	m  = complex(real = cos(id*pi/n), imaginary = sin(id*pi/n));
	return(m);
}


### sum( (zeta(n) - 1) / n )
# Maths 505: A Beautiful Riemann Zeta Sum
# https://www.youtube.com/watch?v=MJn4Fwv3_tU

n = 100
sum(sapply(seq(2, n), \(n) (pracma::zeta(n) - 1) / n))
1 - constEuler


###
n = 200
sum(sapply(seq(2, n), \(n) pracma::zeta(n) / n))
log(n)



### sum( zeta(2*n) * x^(2*n) )
# Michael Penn: Some identities involving the Riemann-Zeta function.
# https://www.youtube.com/watch?v=2W2Ghi9idxM

n = 100
sum(sapply(seq(n), \(n) pracma::zeta(2*n) / 4^n))
1/2

#
sum(sapply(seq(n), \(n) pracma::zeta(2*n) / 2^n))
(1 - pi/sqrt(2) / tan(pi/sqrt(2))) / 2

#
x = 1/sqrt(3)
sum(sapply(seq(n), \(n) pracma::zeta(2*n) * x^(2*n)))
(1 - pi*x/tan(pi*x)) / 2


### sum( (zeta(2*n) - 1) * x^(2*n) )
sum(sapply(seq(n), \(n) pracma::zeta(2*n) - 1))
3/4

#
x = 1/sqrt(3)
sum(sapply(seq(n), \(n) (pracma::zeta(2*n) - 1) * x^(2*n)) )
(3 - pi*x/tan(pi*x)) / 2 - 1/(1 - x^2)


### sum( zeta(4*n) - 1 )
# Michael Penn: More identities involving the Riemann-Zeta function!
# https://www.youtube.com/watch?v=bRdGQKwusiE

n = 500
sum(sapply(seq(n), \(n) zeta(4*n) - 1))
7/8 - pi/4 * (exp(2*pi) + 1) / (exp(2*pi) - 1)

### sum( zeta(6*n) - 1 )
n = 500
sum(sapply(seq(n), \(n) zeta(6*n) - 1))
m3 = cos(2*pi/3) + 1i*c(1,-1)*sin(2*pi/3)
1 - 1/3 + 1/4 - sum(pi*m3/tan(pi*m3)) / 6

### sum( zeta(12*n) - 1 )
n = 500
m12 = unity.neg(12);
sum(sapply(seq(n), \(n) zeta(12*n) - 1))
sum(pi*m12/tan(pi*m12)) / 24 - 1;

### sum( zeta(12*n + 6) - 1 )
n = 300
m3  = cos(2*pi/3) + 1i*c(1,-1)*sin(2*pi/3);
m12 = unity.neg(12);
sum(sapply(seq(0, n), \(n) zeta(12*n + 6) - 1))
2 - 1/3 + 1/4 - sum(pi*m3/tan(pi*m3)) / 6 - sum(pi*m12/tan(pi*m12)) / 24;

### sum( zeta(12*n + 4) + zeta(12*n + 8) - 2 )
n = 300
m12 = unity.neg(12);
sum(sapply(seq(0, n), \(n) zeta(12*n + 4) + zeta(12*n + 8) - 2))
1 + 7/8 - pi/4 * (exp(2*pi) + 1) / (exp(2*pi) - 1) +
	- sum(pi*m12/tan(pi*m12)) / 24;


### sum( zeta(12*n + 4) - 1 )
n = 200
m3  = cos(2*pi/3) + 1i*c(1,-1)*sin(2*pi/3);
m12 = unity.neg(12);
#
sum(sapply(seq(0, n), \(n) zeta(12*n + 4) - 1));
# TODO: simplify
(2*(1 + 7/8 - pi/4 * (exp(2*pi) + 1) / (exp(2*pi) - 1) +
	- sum(pi*m12/tan(pi*m12)) / 24) +
	- m3[1]*(sum((3 - pi*m3[1]*1i/tan(pi*m3[1]*1i)) / 2 - 1/(1 + m3[1]^2)) +
	+ sum((3 - pi*m3[1]/tan(pi*m3[1])) / 2 - 1/(1 - m3[1]^2)) +
	- 2*(sum(pi*m12/tan(pi*m12)) / 24 - 1)) ) / (2 - 2*m3[2]);

### sum( zeta(12*n + 8) - 1 )
sum(sapply(seq(0, n), \(n) zeta(12*n + 8) - 1));
# TODO: simplify
(2*(1 + 7/8 - pi/4 * (exp(2*pi) + 1) / (exp(2*pi) - 1) +
	- sum(pi*m12/tan(pi*m12)) / 24) +
	- m3[2]*(sum((3 - pi*m3[1]*1i/tan(pi*m3[1]*1i)) / 2 - 1/(1 + m3[1]^2)) +
	+ sum((3 - pi*m3[1]/tan(pi*m3[1])) / 2 - 1/(1 - m3[1]^2)) +
	- 2*(sum(pi*m12/tan(pi*m12)) / 24 - 1)) ) / (2 - 2*m3[1]);

# Note:
# sum(zeta(6*n + 2) - 1) = sum(zeta(12*n + 2) + zeta(12*n +  8) - 2)
# sum(zeta(6*n + 4) - 1) = sum(zeta(12*n + 4) + zeta(12*n + 10) - 2)

# Derivation:
sum(sapply(seq(0, n), \(n) {
	m3[2]*(-zeta(12*n + 2) + zeta(12*n + 8)) +
	+ m3[1]*(zeta(12*n + 4) - zeta(12*n + 10))
}))
sum((3 - pi*m3[1]*1i/tan(pi*m3[1]*1i)) / 2 - 1/(1 + m3[1]^2)) +
	- 2*(sum(pi*m12/tan(pi*m12)) / 24 - 1) +
	+ 1 - 1/3 + 1/4 - sum(pi*m3/tan(pi*m3)) / 6;
# for: zeta(12*n + 4) & zeta(12*n + 8)
sum(sapply(seq(0, n), \(n) {
	m3[2]*zeta(12*n + 8) + m3[1]*zeta(12*n + 4) + 1;
})) * 2
sum((3 - pi*m3[1]*1i/tan(pi*m3[1]*1i)) / 2 - 1/(1 + m3[1]^2)) +
	+ sum((3 - pi*m3[1]/tan(pi*m3[1])) / 2 - 1/(1 - m3[1]^2)) +
	- 2*(sum(pi*m12/tan(pi*m12)) / 24 - 1);

# for: zeta(12*n + 2) & zeta(12*n + 10)
sum(sapply(seq(0, n), \(n) {
	m3[2]*zeta(12*n + 2) + m3[1]*zeta(12*n + 10) + 1;
})) * -2
sum((3 - pi*m3[1]*1i/tan(pi*m3[1]*1i)) / 2 - 1/(1 + m3[1]^2)) +
	- sum((3 - pi*m3[1]/tan(pi*m3[1])) / 2 - 1/(1 - m3[1]^2)) +
	+ 2*(2 - 1/3 + 1/4 - sum(pi*m3/tan(pi*m3)) / 6 - sum(pi*m12/tan(pi*m12)) / 24);


##################

###
# Michael Penn: Another Riemann-Zeta function identity.
# https://www.youtube.com/watch?v=AmWvTN6cPz8

n = 500
sum(sapply(seq(n), \(n) pracma::zeta(2*n) / ((n+1)*(2*n+1))) )
1/2
#
n = 2000
sum(sapply(seq(n), \(n) zeta(2*n) / ((n+1)*(2*n+1))) )
1/2


###
# Michael Penn: More Riemann Zeta function identities!!
# https://www.youtube.com/watch?v=JwxgwXUruRM

n = 500
sum(sapply(seq(n), \(n) zeta(2*n) - zeta(2*n+1)) )
1/2

###
n = 500
sum(sapply(seq(n), \(n) zeta(2*n+1) - 1) )
1/4

###
n = 1000
sum(sapply(seq(n), \(n) zeta(n+1) - 1) )
1


#################

###
# Maths 505: 3 RIDICULOUSLY AWESOME infinite zeta series!!!
# https://www.youtube.com/watch?v=Pl-GUHN4DyI

# sum[n]( zeta(n+1)*x^n )
# = sum[k]( 1/(k-x) - 1/k )
# = - digamma(1-x) - constEuler;

n = 500

###
x = 2
sum(sapply(seq(n), \(n) zeta(n+1) / x^n))
- digamma(1 - 1/x) - constEuler;

###
x = sqrt(5)
sum(sapply(seq(n), \(n) zeta(n+1) / x^n))
- digamma(1 - 1/x) - constEuler;


### sum( zeta(n) / x^n )
# Note: from n = 2!
sum(sapply(seq(2, 100), \(n) zeta(n) / 2^n))
log(2)

###
x = 3
sum(sapply(seq(2, 100), \(n) zeta(n) / x^n))
- (digamma(1 - 1/x) + constEuler) / x

### sum( zeta(2*n+1) / x^(2*n+1) )
# =>
x = sqrt(5)
sum(sapply(seq(100), \(n) zeta(2*n+1) / x^(2*n+1)))
- (digamma(1 - 1/x) + constEuler) / x - (1 - pi/x / tan(pi/x)) / 2;
- (digamma(1 - 1/x) + constEuler - pi/2/tan(pi/x)) / x - 1/2;
- (digamma(1 + 1/x) + digamma(1 - 1/x) + 2*constEuler) / (2*x);
- (digamma(1/x) + digamma(- 1/x) + 2*constEuler) / (2*x);


### sum( (-1)^n * zeta(n+1) / x^n )
x = 2
sum(sapply(seq(100), \(n) - (-1)^n * zeta(n+1) / x^n))
digamma(1 + 1/x) + constEuler

###
x = sqrt(3)
sum(sapply(seq(100), \(n) - (-1)^n * zeta(n+1) / x^n))
digamma(1 + 1/x) + constEuler


### sum( (-1)^n * zeta(n) / x^n )
# from n = 2;
x = 3
sum(sapply(seq(2, 100), \(n) (-1)^n * zeta(n) / x^n))
(digamma(1 + 1/x) + constEuler) / x


### sum( n * zeta(n+1) / x^n )
# = sum( 1 / (k - 1/x)^2 ) / x;

x = 2
sum(sapply(seq(n), \(n) n * zeta(n+1) / x^n))
pracma::psi(1, 1 - 1/x) / x

###
x = sqrt(5)
sum(sapply(seq(n), \(n) n * zeta(n+1) / x^n))
pracma::psi(1, 1 - 1/x) / x


#################

### Varia

###
# only Integers!
x = 7
digamma(1/x) + digamma(-1/x)
id = seq(floor(x/2));
cs = cos(2*pi*id/x); sn = sin(pi*id/x);
x - 2*constEuler - 2*log(2*x) + 4*sum(cs * log(sn))

###
x = sqrt(7)
digamma(1/x) - digamma(-1/x)
- x - pi/tan(pi/x)

###
digamma(-1/x)
digamma(1/x) + pi/tan(pi/x) + x


######################
######################

### sum( (-1)^n * zeta(n) * x^n / n )
# Maths 505: THE MOST BEAUTIFUL RESULT IN ALL OF CALCULUS!!!
# https://www.youtube.com/watch?v=9BTiRweEcT0


###
id = seq(2, 1000)
sum(zeta(id) / (id * 2^id))
log(pi / exp(Euler)) / 2

###
x = -1/3
id = seq(2, 1000)
log(gamma(x + 1))
sum((-1)^id * zeta(id) * x^id / id) - Euler * x;

###
x = 1/sqrt(3); # x = 1/pi;
id = seq(1, 1000)
log(pi*x/sin(pi*x))
sum(zeta(2*id) * x^(2*id) / id);

