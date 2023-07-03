

### S5 Ht Theory

# - the start to develop a Theory;


##################

### Helper

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

x = sqrt(c(3,5,7,11,13));
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

E2a = as.E2ps(x);
E3a = as.E3ps(x, p=1);
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


### Prod( x - x1*x2*x3*x4 )
# [symmetric]
E4p5 - E4*E4p4 + E5*E3*E4p3 - E5^2*E2*E4p2 + E5^3*S*E4 - 5*E5^4


### Ea-based:

### Ht Prod( x - x1*x2*x3 )
E3a5 - E3a*E3a4 + (E5*S + E1221a)*E3a3 - E5*(E4 + E121a)*E3a2 + E5^2*E2a*E3a - 5*E5^3

### Ht Prod( x - x1*x2 )
E2a5 - E2a*E2a4 + (E4 + E121a)*E2a3 - (E5*S + E1221a)*E2a2 + E5*E3a*E2a - 5*E5^2


