########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S6: 2 x C3-Hetero-Symmetric
### Basic Types
###
### draft v.0.1a


####################

### Helper Functions


source("Polynomials.Helper.R")


#####################
#####################

# for other Types: see Poly.System.S6.C3.R;

### System:
# x1 + x2 + x3 = Sx
# y1 + y2 + y3 = Sy
# x1*x2 + x2*x3 + x3*x1 = E2x
# y1*y2 + y2*y3 + y3*y1 = E2y
# x1*y1 + x2*y2 + x3*y3 = R5
# px * py = E6

### Note:
# - IF {(x1, x2, x3), (y1, y2, y3)} is a solution,
#   and sigma is a given permutation, then:
#   {sigma(x1, x2, x3), sigma(y1, y2, y3)} is also a solution;


### Solver

solve.S6.2C3 = function(R, debug = TRUE) {
	coeff = coef.S6.2C3(R);
	E6 = R[6];
	px = roots(coeff);
	py = E6 / px;
	if(debug) print(px);
	sx = R[1]; sy = R[2];
	E2x = R[3]; E2y = R[4];
	x1 = sapply(px, function(px) {
		roots(c(1, - sx, E2x, - px));
	})
	# TODO
	return(x1);
}
coef.S6.2C3 = function(R) {
	sx = R[1]; sy = R[2];
	E2x = R[3]; E2y = R[4];
	E2a = R[5]; E6 = R[6];
	c((sy^2 - 3*E2y)^3,
	sy^2*E2y^2*sx^3 - 4*E2y^3*sx^3 + sy^4*E2y*sx*E2x - 9*sy^2*E2y^2*sx*E2x + 18*E2y^3*sx*E2x +
		- 2*sy^3*E2y*sx^2*E2a + 9*sy*E2y^2*sx^2*E2a - 2*sy^5*E2x*E2a + 15*sy^3*E2y*E2x*E2a - 27*sy*E2y^2*E2x*E2a +
		+ 2*sy^4*sx*E2a^2 - 9*sy^2*E2y*sx*E2a^2 - 2*sy^3*E2a^3 + 9*sy*E2y*E2a^3,
	E2y^3*sx^2*E2x^2 + sy^2*E2y^2*E2x^3 - 4*E2y^3*E2x^3 - sy*E2y^2*sx^3*E2x*E2a - sy^3*E2y*sx*E2x^2*E2a +
		+ 3*sy*E2y^2*sx*E2x^2*E2a + E2y^2*sx^4*E2a^2 + 3*sy^2*E2y*sx^2*E2x*E2a^2 - 6*E2y^2*sx^2*E2x*E2a^2 +
		+ sy^4*E2x^2*E2a^2 - 6*sy^2*E2y*E2x^2*E2a^2 + 9*E2y^2*E2x^2*E2a^2 - 2*sy*E2y*sx^3*E2a^3 +
		- 2*sy^3*sx*E2x*E2a^3 + 5*sy*E2y*sx*E2x*E2a^3 + sy^2*sx^2*E2a^4 + 2*E2y*sx^2*E2a^4 + 2*sy^2*E2x*E2a^4 +
		- 6*E2y*E2x*E2a^4 - 2*sy*sx*E2a^5 + E2a^6 - 2*sy^3*sx^3*E6 + 9*sy*E2y*sx^3*E6 + 9*sy^3*sx*E2x*E6 +
		- 27*sy*E2y*sx*E2x*E6 - 27*E2y*sx^2*E2a*E6 - 27*sy^2*E2x*E2a*E6 + 81*E2y*E2x*E2a*E6 + 27*sy*sx*E2a^2*E6 +
		- 27*E2a^3*E6,
	sy*E2y*sx^4*E2x*E6 + sy^3*sx^2*E2x^2*E6 - 9*sy*E2y*sx^2*E2x^2*E6 - 4*sy^3*E2x^3*E6 +
		+ 18*sy*E2y*E2x^3*E6 - 2*E2y*sx^5*E2a*E6 - 2*sy^2*sx^3*E2x*E2a*E6 + 15*E2y*sx^3*E2x*E2a*E6
		+ 9*sy^2*sx*E2x^2*E2a*E6 - 27*E2y*sx*E2x^2*E2a*E6 + 2*sy*sx^4*E2a^2*E6 - 9*sy*sx^2*E2x*E2a^2*E6 +
		- 2*sx^3*E2a^3*E6 + 9*sx*E2x*E2a^3*E6,
	(sx^2 - 3*E2x)^3*E6^2);
}

### Test
x = sqrt(c(2,3,5,7,11,13))
y = x[4:6]; x = x[1:3];
sx  = sum(x); sy = sum(y);
E2x = x[1]*x[2] + x[1]*x[3] + x[2]*x[3];
E2y = y[1]*y[2] + y[1]*y[3] + y[2]*y[3];
px  = prod(x); py = prod(y); E6 = px * py;
E2a = sum(x*y); # R5
E2b = sum(x * y[c(2,3,1)]);
E2c = sum(x * y[c(3,1,2)]);
E2abc = E2a*(E2b + E2c) + E2b*E2c;
# [intermediary]
A2a = sum(x^2 * x[c(2,3,1)]);
A2b = sum(x^2 * x[c(3,1,2)]);
B2a = sum(y^2 * y[c(2,3,1)]);
B2b = sum(y^2 * y[c(3,1,2)]);

# Test
x^3 - sx*x^2 + E2x*x - px # = 0
y^3 - sy*y^2 + E2y*y - py # = 0
sum(coef.S6.2C3(c(sx, sy, E2x, E2y, E2a, E6)) * px^seq(4, 0)) # = 0

### Derivation:

# E2b + E2c:
E2a + E2b + E2c - sx*sy # = 0

# E2abc = E2 for (E2a, E2b, E2c)
E2abc + 3*E2x*E2y - E2y*sx^2 - E2x*sy^2 # = 0
# =>
E2b*E2c + E2a*(sx*sy - E2a) + 3*E2x*E2y - E2y*sx^2 - E2x*sy^2 # = 0
# => known E2a; solve E2b, E2c;

# E3 for (E2a, E2b, E2c)
E2a*E2b*E2c - py*(sx^3 - 3*E2x*sx) - px*(sy^3 - 3*E2y*sy) - 18*px*py +
	- 2*A2a*B2a + (E2x*sx - 3*px)*B2a + (E2y*sy - 3*py)*A2a +
	- E2x*sx*E2y*sy + 3*E2x*py*sx + 3*E2y*px*sy # = 0
# =>
E2a*(E2y*sx^2 + E2x*sy^2 - 3*E2x*E2y - E2a*(sx*sy - E2a)) +
	- py*(sx^3 - 3*E2x*sx) - px*(sy^3 - 3*E2y*sy) - 18*px*py +
	- 2*A2a*B2a + (E2x*sx - 3*px)*B2a + (E2y*sy - 3*py)*A2a +
	- E2x*sx*E2y*sy + 3*E2x*py*sx + 3*E2y*px*sy # = 0
# with:
A2a^2 - (E2x*sx - 3*px)*A2a - (6*E2x*px*sx - E2x^3 - 9*px^2 - px*sx^3) # = 0
B2a^2 - (E2y*sy - 3*py)*B2a - (6*E2y*py*sy - E2y^3 - 9*py^2 - py*sy^3) # = 0

pE3 = as.pm("...") # see formula for E3 above
pA2 = as.pm("A2a^2 - (E2x*sx - 3*px)*A2a - (6*E2x*px*sx - E2x^3 - 9*px^2 - px*sx^3)")
pB2 = as.pm("B2a^2 - (E2y*sy - 3*py)*B2a - (6*E2y*py*sy - E2y^3 - 9*py^2 - py*sy^3)")

pR = solve.pm(pA2, pE3, by = "A2a")
pR = solve.pm(pB2, pR$Rez, by = "B2a")
pR = pR$Rez;
pR = replace.fr.pm(pR, as.pm("E6"), as.pm("px"), "py")
as.coef.pm(pR, "px")


