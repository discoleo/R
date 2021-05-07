########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S2:
### Mixed Type: Composite
###
### draft v.0.1b


### Heterogeneous Symmetric
### Polynomial Systems: 2 Variables
### Mixed/Composite: Hetero + Symmetric

### Example:
# x^p + b*y = Ru # where Ru = unknown
# y^p + b*x = Ru
# x*y = R2

### Variants:
# x^n + y^n = R2
# x*y*(x+y) = R2
# or
# another Hetero-Symmetric system;


######################
######################

###############
### History ###
###############


### draft v.0.1b:
# - Order 3 variant:
#   x*y*(x+y) = R;
### draft v.0.1a:
# - solved: double Ht-Symmetric Composite;
# - example P[6]:
#   3 - 2*x + x^2 + 4*x^3 - x^4 + x^6 = 0;


######################
######################

################
### Order: 3 ###
################

# x^3 + b2*y^2 + b1*y = Ru1
# y^3 + b2*x^2 + b1*x = Ru1
# x*y*(x+y) = R

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^2 - x*y - b2*S - b1 # = 0
# x*y = S^2 - b2*S - b1;

### =>
S^3 - b2*S^2 - b1*S - R # = 0


### Solver:
solve.HtComposite.S2P3 = function(b, debug=TRUE) {
	coeff = c(1, -b[2], -b[1], -R[1]);
	S = roots(coeff);
	if(debug) print(S);
	xy = S^2 - b[2]*S - b[1];
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	sol = cbind(x=x, y=y);
	sol = rbind(sol, sol[,2:1]);
	return(sol);
}

### Examples:

R = -1
b = c(1, -2)
sol = solve.HtComposite.S2P3(b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^3 + b[2]*y^2 + b[1]*y, y^3 + b[2]*x^2 + b[1]*x)
x*y*(x+y) # - R[1]

###
round0.p(poly.calc(x))
err = 1 + 6*x - x^2 - 4*x^3 - 2*x^4 + 2*x^5 + x^6
round0(err)


######################

######################
### Full Composite ###
######################

####################
### Order: 3 & 4 ###
####################

# x^3 + b1*y = Ru1
# y^3 + b1*x = Ru1
# x^4 + b2*y = Ru2
# y^4 + b2*x = Ru2

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^2 - x*y - b1 # = 0
# x*y = S^2 - b1;

### Diff Eqs 3, 4 =>
S*(S^2 - 2*x*y) - b2 # = 0

### =>
S*(S^2 - 2*b1) + b2 # = 0
S^3 - 2*b1*S + b2 # = 0

### Solver:
solve.HtComposite2.S2P34 = function(b, debug=TRUE) {
	coeff = c(1, 0, -2*b[1], b[2]);
	S = roots(coeff);
	if(debug) print(S);
	xy = S^2 - b[1];
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	sol = cbind(x=x, y=y);
	sol = rbind(sol, sol[,2:1]);
	return(sol);
}

### Examples:

b = c(1, -2)
sol = solve.HtComposite2.S2P34(b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^3 + b[1]*y, y^3 + b[1]*x)
cbind(x^4 + b[2]*y, y^4 + b[2]*x)

###
round0.p(poly.calc(x))
err = 3 - 2*x + x^2 + 4*x^3 - x^4 + x^6
round0(err)

