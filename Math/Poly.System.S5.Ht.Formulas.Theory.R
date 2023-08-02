

### S5 Ht Theory

# - the start to develop a Theory;


##################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


shift.f = function(x, shift = 1) {
	if(shift == 0) return(x);
	c(x[-seq(shift)], x[seq(shift)]);
}
# symmetric
as.E2ps = function(x, pow = 1, shift = 1) {
	if(pow != 1) x = x^pow;
	sum(x * shift.f(x, shift=shift));
}
# symmetric
as.E3ps = function(x, pow = 1) {
	if(pow != 1) x = x^pow;
	sum(x * c(x[-1], x[1]) * c(x[-c(1,2)], x[c(1,2)]));
}
as.E3p = function(x, pow=c(1,2,1)) {
	y = x^pow[2]; z = x^pow[3]; x = x^pow[1];
	sum(x * c(y[-1], y[1]) * c(z[-c(1,2)], z[c(1,2)]));
}
as.E4 = function(x, pow = 1) {
	x = x^pow;
	s1 = x[1] + x[3]; s2 = x[2] + x[4];
	p1 = x[1] * x[3]; p2 = x[2] * x[4];
	p1*p2 + x[5]*(s1*p2 + s2*p1);
}
as.E4p = function(x, pow = c(1,2,2,1)) {
	x2 = x^pow[2]; x3 = x^pow[3];
	x4 = x^pow[4]; x1 = x^pow[1];
	sum(x1 * shift.f(x2,1) * shift.f(x3,2) * shift.f(x4,3));
}


#################

### Test Values

x = (c(3,5,7,11,13))^(1/3);
x[4] = - x[4];

s1 = x[1] + x[3]; s2 = x[2] + x[4];
p1 = x[1] * x[3]; p2 = x[2] * x[4];
S  = s1 + s2 + x[5]; E2p = s1 * s2 + p1 + p2;
E2 = x[5] * (s1 + s2) + E2p;
E3 = E2p * x[5] + s1*p2 + s2*p1;
E4 = p1 * p2 + x[5]*(s1*p2 + s2*p1);
E5 = p1 * p2 * x[5];

# Test:
x^5 - S*x^4 + E2*x^3 - E3*x^2 + E4*x - E5 # == 0

S2  = sum(x^2);
E2a = as.E2ps(x);
E3a = as.E3ps(x, p=1);
E2b = E2 - E2a;
E22a = as.E2ps(x, pow=2);
E22b = as.E2ps(x, pow=2, shift=2);
#
E121a  = as.E3p(x, pow=c(1,2,1));
E1221a = as.E4p(x, pow=c(1,2,2,1));

# Powers:
E4p5 = as.E4(x, 5); E4p4 = as.E4(x, 4);
E4p3 = as.E4(x, 3); E4p2 = as.E4(x, 2);
#
E3a2 = as.E3ps(x, p=2); E3a3 = as.E3ps(x, p=3);
E3a4 = as.E3ps(x, p=4); E3a5 = as.E3ps(x, p=5);
#
E2a2 = as.E2ps(x, p=2); E2a3 = as.E2ps(x, p=3);
E2a4 = as.E2ps(x, p=4); E2a5 = as.E2ps(x, p=5);


#######################

###
# D. Aharonov, U. Elias. More on the identity of Chaundy and Bullard. (2014)
# https://doi.org/10.1016/j.jmaa.2014.04.025

###
1 / E5
sum( x / (S*E5) )
sum( x * shift.f(x, 1) / (E2a*E5) )
sum( x * shift.f(x, 2) / (E2b*E5) )

# TODO:
# - find way to combine all components: S, E2a, E2b, E3, E4, E5;


#####################

### Prod( x - x1*x2*x3*x4 )
# [symmetric]
E4p5 - E4*E4p4 + E5*E3*E4p3 - E5^2*E2*E4p2 + E5^3*S*E4 - 5*E5^4


### Ea-based:

### Ht Prod( x - x1*x2*x3 )
E3a5 - E3a*E3a4 + (E5*S + E1221a)*E3a3 - E5*(E4 + E121a)*E3a2 + E5^2*E2a*E3a - 5*E5^3

### Ht Prod( x - x1*x2 )
E2a5 - E2a*E2a4 + (E4 + E121a)*E2a3 - (E5*S + E1221a)*E2a2 + E5*E3a*E2a - 5*E5^2

### Additional Eqs:
E2a2 + 2*E121a + 2*E4 - E2a^2 # = 0
E3a2 + 2*E1221a + 2*E5*S  - E3a^2 # = 0

# TODO: remaining;


### Ht Prod( x - x1 - x2 )
r = x + shift.f(x, 1)
poly.calc(r)

# Coeff of x^4:
- 2*S
# Coeff of x^3:
sum(x^2, 4*E2, - E2a)
# Coeff of x^2:
sum(6*E3, -2*E3a) + S*sum(x^2) - sum(x^3) + sum(x^2 * (shift.f(x,2) + shift.f(x,3)))
7*E3 - E3a - E2a*S + 2*S*(S^2 - 2*E2) - 2*sum(x^3)
# Coeff of x:
E121a - E3a*S - E2a*E2 + E3*S + E2^2

# TODO: Coeff of b0;

### Powers of x:
sum(- r^2, 2*x^2, 2*E2a)
sum(- r^3, 2*x^3, 3*E2a*S, -3*E3, -3*E3a)
sum(- r^4, 2*x^4, 2*E2a^2 + 4*(E2a*(S^2 - E2) - E121a - E3a*S - E3*S + 2*E4))

# TODO: r^5;


### Ht Prod( x - x1 - x3 )
r = x + shift.f(x, 2)
poly.calc(r)

# TODO


### Helper & Derivation:
S2^2 - sum(x^4) - 2*E22a - 2*E22b # = 0
E2a^2 - E22a - 2*E121a - 2*E4 # = 0
# E2101a + E1201a
E2a*E2b - E3a*S + E121a + 2*E4 - E4 +
	- sum(x^2 * (shift.f(x,1) * shift.f(x,3) + shift.f(x,4) * shift.f(x,2))) # = 0

# E31a + E13a
sum(x^3 * (shift.f(x,1) + shift.f(x,4)))
sum(E2a*x^2, - E3a*S, E121a, 2*E4, - x * shift.f(x,1) * shift.f(x^2,3))
# =>
sum(E2a*x^2, - E3a*S, E121a, 2*E4) + (E2a*E2b - E3a*S + E121a + E4) +
	- sum(x * shift.f(x,1) * shift.f(x^2,3)) +
	- sum(x^2 * (shift.f(x,1) * shift.f(x,3) + shift.f(x,4) * shift.f(x,2))) # = 0
sum(E2a*x^2) + E2a*E2b + 2*E121a - E3a*S - E3*S + 5*E4
E2a*E2b - 2*E2*E2a + E2a*S^2 + 2*E121a - E3a*S - E3*S + 5*E4

# =>
r = x + shift.f(x, 1)
sum(- r^4, 2*x^4, 4 * sum(x^3 * (shift.f(x,1) + shift.f(x,4))),
	6 * sum(x^2 * shift.f(x^2,1)) )
sum(- r^4, 2*x^4, 4 * E2a*x^2, 4*(E2a*E2b + 2*E121a - E3a*S - E3*S + 5*E4),
	6 * sum(x^2 * shift.f(x^2,1)) )
sum(- r^4, 2*x^4, 4*E2a*x^2, 6*E2a^2 + 4*(E2a*E2b - E121a - E3a*S - E3*S + 2*E4))
sum(- r^4, 2*x^4, 2*E2a^2 + 4*(E2a*(S^2 - E2) - E121a - E3a*S - E3*S + 2*E4))


### Derivation: Coeff of x
poly.calc(r)
2*sum((x^2*shift.f(x,1) + x*shift.f(x^2,1)) * shift.f(x,3)) + # E2101a + E1201a
	+ sum(x^2*shift.f(x,1) * shift.f(x,2),
		x*shift.f(x,1) * shift.f(x^2,2)) + # E211a + E112a
	+ 3*sum(x^2*shift.f(x,2) * shift.f(x,3)) + # E2011a
	+ 5*E4 + E121a + sum(x^2 * shift.f(x^2, 2))
2*sum((x^2*shift.f(x,1) + x*shift.f(x^2,1)) * shift.f(x,3)) + # E2101a + E1201a
	+ 3*sum(x^2*shift.f(x,2) * shift.f(x,3)) + # E2011a
	+ 3*E4 + E3a*S + sum(x^2 * shift.f(x^2, 2))
2*E3*S - E3a*S - E4 +
	+ sum(x^2 * shift.f(x,2) * shift.f(x,3)) + # E2011a
	+ sum(x^2 * shift.f(x^2, 2))
2*E3*S - E3a*S - E4 - (E2a*E2b + E121a - E3*S + 3*E4) +
	+ sum(x^2 * shift.f(x^2, 2))
E121a - E3a*S - E2a*E2 + E3*S + E2^2

### E2011a
sum(x^2 * shift.f(x,2) * shift.f(x,3))
- (E2a*E2b + E121a - E3*S + 3*E4)

### E202a
sum(x^2 * shift.f(x^2, 2))
S2^2/2 - sum(x^4)/2 + (2*E121a + 2*E4 - E2a^2)
2*E121a + 4*E4 - 2*E3*S + E2^2 - E2a^2


