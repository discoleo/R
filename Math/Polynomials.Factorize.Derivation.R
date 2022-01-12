########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Factorize: Derivations
###
### draft v.0.1b-fix3


### Factorize Multi-Variable Polynomials
### Various Derivations

### A. Univariate Polynomials


######################

### Helper Functions

source("Polynomials.Helper.R")
# - is automatically loaded in: Polynomials.Helper.R;
# source("Polynomials.Helper.Factorize.R")


#######################
#######################

### A. Univariate Polynomials
### P[3] * P[3]

### Factors: (x^3 + ... + 1)
# (x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1)
p = toPoly.pm("(x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1)");

b = c(2, 4); c = c(3, -3);
b1 = b[1]; b2 = b[2]; c1 = c[1]; c2 = c[2];
P = function(x) eval.pm(p, list(x=x, b1=b1, b2=b2, c1=c1, c2=c2));

# TODO: recheck calculations;

#########
### Mod 7

### P(3) (mod 7)
(9*(3*b2 + b1)*(3*c2 + c1) - P(3)) %% 7
((3*b2 + b1)*(3*c2 + c1) + 3*P(3)) %% 7
(b1*c1 + 2*b2*c2 + 3*(b1*c2 + b2*c1) + 3*P(3)) %% 7

### P(5) (mod 7)
(25*(5*b2 + b1)*(5*c2 + c1) - P(5)) %% 7
((5*b2 + b1)*(5*c2 + c1) - 2*P(5)) %% 7
(b1*c1 + 4*b2*c2 + 5*(b1*c2 + b2*c1)- 2*P(5)) %% 7

### P(6) (mod 7)
((b2 - b1)*(c2 - c1) - P(6)) %% 7

### Aux: P(1) (mod 7)
((b2 + b1 + 2)*(c2 + c1 + 2) - P(1)) %% 7
### Aux: P(2) (mod 7)
(4*(2*b2 + b1 + 1)*(2*c2 + c1 + 1) - P(2)) %% 7
((2*b2 + b1 + 1)*(2*c2 + c1 + 1) - 2*P(2)) %% 7

# System =>
(4*b1*c1 + 5*b2*c2 + 3*P(3) - 3*P(6)) %% 7
(6*b1*c1 + 2*b2*c2 - 2*P(5) - 5*P(6)) %% 7
# =>
(b1*c1 + 3*b2*c2 - P(3) + P(6)) %% 7
(b1*c1 - 2*b2*c2 + 2*P(5) + 5*P(6)) %% 7
# Eqs:
(b1*c1 + P(3) - 3*P(5) + 2*P(6)) %% 7
(b2*c2 - 3*P(3) + P(5) + 2*P(6)) %% 7
((b1*c2 + b2*c1) - 2*P(3) - 2*P(5) + 5*P(6)) %% 7
# Aux =>
(b1*c1 + b2*c2 + (b1*c2 + b2*c1) + 2*(b1+b2+c1+c2) + 4 - P(1)) %% 7
(2*(b1+b2+c1+c2) + 4 - P(1) + 4*P(3) + 4*P(5) - 2*P(6)) %% 7
((b1+b2) + (c1+c2) + 2 - 4*P(1) + 2*P(3) + 2*P(5) - P(6)) %% 7
# =>
((b1+b2)*(c1+c2) - 4*P(3) - 4*P(5) + 2*P(6)) %% 7
### Eq:
x = (b1 + b2);
(x^2 + (2 - 4*P(1) + 2*P(3) + 2*P(5) - P(6))*x + 4*P(3) + 4*P(5) - 2*P(6)) %% 7

b1.f = function(P) (2 - 4*P(1) + 2*P(3) + 2*P(5) - P(6)) %% 7;
b0.f = function(P) (4*P(3) + 4*P(5) - 2*P(6)) %% 7;

### Example P1:
b1.f(P); b0.f(P);
# x^2 + x + 0 = 0 =>
# x =  0: c1 + c2 = 0 OR
# x = -1: b1 + b2 = -1;
# =>
b1 + b2 + 1 # = 0 (Mod 7)
c1 + c2 # = 0
b1*c1 + 1 # = 0
b2*c2 - 2 # = 0
# =>
(1 + b1)*(c1) - 2 # = 0 (Mod 7)
# => (c1 + (b1*c1) - 2) %% 7 # = 0, where b1*c1 = -1;
c(c1 - 3, c2 + 3) %% 7 # = 0
c(b1 - 2, b2 - 4) %% 7 # = 0
# b2*c2 - 2 == (-14);


### Example P2: NO Solution
P2 = function(x, b1=2) x^6 + b1*x + 1;
b1.f(P2); b0.f(P2);
# x^2 + 2*x + 4 = 0
sapply(seq(7), function(x) (x^2 + 2*x + 4) %% 7)
x = c(1, 4)
((x - 1)*(x - 4)) %% 7
# =>
b1 + b2 - 1 # = 0
c1 + c2 - 4 # = 0
b1*c1 + 0 # = 0
b2*c2 + 2 # = 0
# =>
(b1 - 1)*(c1 - 4) + 2 # = 0
b1*c1 - 4*b1 - c1 + 6 # = 0
4*b1 + c1 - 6 # = 0
# =>
4*b1^2 - 6*b1 # = 0
b1*(b1 + 2) # = 0
#
b1 = c(0, 5);
b2 = c(1, 3);
c1 = c(6, 0);
c2 = c(5, 4);

### Aux:
((2*b2 + b1 + 1)*(2*c2 + c1 + 1) - 2*P2(2)) %% 7 == 0 # FALSE !
(b1*c1 + 4*b2*c2 + 2*(b1*c2+b2*c1) + (b1+c1) + 2*(b2+c2) + 1 - 2*P2(2)) %% 7
# => NOT Factorisable!

solve.ModP2 = function(b, mod) {
	b0 = b[1]; b1 = b[2]; b2 = b[3];
	x.f = function(x) (b2*x^2 + b1*x + b0) %% mod;
	if(b0 == 0) {
		if(b2 == 1) {
			x = c(0, (mod - b1));
		} else if(b2 == 0) {
			x = 0;
		} else {
			b2inv = inv.mod(b2, mod=mod);
			x = c(0, (- b1*b2inv) %% mod);
		}
	} else {
		err = sapply(seq(mod), x.f);
		x   = which(err == 0);
		if(length(x) == 0) return(list(hasSol = FALSE, Mod = mod));
	}
	return(list(hasSol = TRUE, Sol = x, Mod = mod));
}
inv.mod = function(x, mod) {
	if(mod == 1 || mod == 0) stop("Invalid mod!");
	f = function(x) {
		id = which((x * seq(mod - 1) - 1) %% mod == 0);
		if(length(id) == 0) NA else id;
	}
	sapply(x, f);
}
factorize.S1P6.F3F3 = function(p, debug=FALSE) {
	mod = 7;
	p1x = eval.pm(p, 1); p3x = eval.pm(p, 3);
	p5x = eval.pm(p, 5); p6x = eval.pm(p, 6);
	p1 = p1x %% 7; p3 = p3x %% 7;
	p5 = p5x %% 7; p6 = p6x %% 7;
	b1x = (2 - 4*p1 + 2*p3 + 2*p5 - p6) %% 7;
	b0x = (4*p3 + 4*p5 - 2*p6) %% 7;
	#
	sol = solve.ModP2(c(b0x, b1x, 1), mod=mod);
	if( ! sol$hasSol) return(list(isF = FALSE, Mod = mod));
	b12s = sol$Sol;
	c12s = (- b1x - b12s) %% mod;
	b1c1 = (- p3 + 3*p5 - 2*p6) %% mod;
	b2c2 = (3*p3 - p5 - 2*p6) %% mod;
	len  = length(b12s);
	b1c1 = rep(b1c1, len); b2c2 = rep(b2c2, len);
	if(debug) {
		cat(c("\nEq 1: ", c(b1x, b0x)));
		cat(c("\nb1 + b2 = ", b12s));
		cat(c("\nc1 + c2 = ", c12s));
		cat(c("\nb1 * c1 = ", b1c1));
		cat(c("\nb2 * c2 = ", b2c2, "\n"));
	}
	# (b1 - b12s)*(c1 - c12s) - b2c2 = 0
	# c12s*b1^2 - (b1c1 - b2c2 + b12s*c12s)*b1 + b1c1*b12s = 0
	b1y = (- (b1c1 - b2c2 + b12s*c12s)) %% mod;
	b0y = (b1c1*b12s) %% mod;
	len = length(b0y);
	sol = lapply(seq(len), function(id) solve.ModP2(c(b0y[id], b1y[id], c12s[id]), mod=mod));
	# Filter Solutions:
	len = length(sol);
	b1  = unlist(lapply(sol, function(sol) if(sol$hasSol) sol$Sol else NULL));
	if(length(b1) == 0) return(list(isF = FALSE, mod = mod));
	#
	nSols = unlist(lapply(seq(len), function(id) if(sol[[id]]$hasSol) rep(id, length(sol[[id]]$Sol)) else c()));
	# print(nSols); print(b1);
	b12s = b12s[nSols]; c12s = c12s[nSols];
	b1c1 = b1c1[nSols]; b2c2 = b2c2[nSols];
	b2 = (b12s - b1) %% mod;
	# c1 & c2:
	len = length(b1)
	c1 = rep(NA, len); c2 = rep(NA, len);
	for(id in seq(len)) {
		if(b1[id] != 0) {
			c1[id] = (b1c1[id] * inv.mod(b1[id], mod=mod)) %% mod;
			c2[id] = (c12s[id] - c1[id]) %% mod;
		} else {
			c2[id] = (b2c2[id] * inv.mod(b2[id], mod=mod)) %% mod;
			c1[id] = (c12s[id] - c2[id]) %% mod;
		}
	}
	sol = cbind(b1=b1, b2=b2, c1=c1, c2=c2);
	isNA = apply(sol, 1, function(x) any(is.na(x)));
	sol  = sol[ ! isNA, , drop=FALSE];
	if(nrow(sol) == 0) return(list(isF = FALSE, Mod=mod));
	### Aux Eqs:
	if(debug) {
		print(paste0("Starting Aux eqs: sol = ", nrow(sol)));
		print(sol);
	}
	p2 = eval.pm(p, 2) %% mod;
	p4 = eval.pm(p, 4) %% mod;
	b1 = sol[,1]; b2 = sol[,2]; c1 = sol[,3]; c2 = sol[,4];
	# P(2)
	err1 = ((2*b2 + b1 + 1)*(2*c2 + c1 + 1) - 2*p2) %% mod;
	# P(4)
	err2 = ((b2 + 2*b1 + 1)*(c2 + 2*c1 + 1) - 2*p4) %% mod;
	isOK = (err1 == 0) & (err2 == 0);
	sol = sol[isOK, , drop=FALSE];
	if(nrow(sol) == 0) return(list(isF = FALSE, Mod=mod));
	return(list(isF = TRUE, F=sol, Mod=mod));
}

###
p = toPoly.pm("(x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1)");

###
b = c(2, 4); c = c(3, -6);
b1 = b[1]; b2 = b[2]; c1 = c[1]; c2 = c[2];
p2 = replace.pm(p, list(b1=b1, b2=b2, c1=c1, c2=c2))
P = function(x) eval.pm(p2, x);
factorize.S1P6.F3F3(p2)

###
b = c(2, 2); c = c(11, 12);
b1 = b[1]; b2 = b[2]; c1 = c[1]; c2 = c[2];
p2 = replace.pm(p, list(b1=b1, b2=b2, c1=c1, c2=c2))
P = function(x) eval.pm(p2, x);
factorize.S1P6.F3F3(p2)

###
b = c(8, 7); c = c(-11, 8);
b1 = b[1]; b2 = b[2]; c1 = c[1]; c2 = c[2];
p2 = replace.pm(p, list(b1=b1, b2=b2, c1=c1, c2=c2))
P = function(x) eval.pm(p2, x);
factorize.S1P6.F3F3(p2)

### False-Positives:
# What is the meaning?
factorize.S1P6.F3F3(toPoly.pm("x^6 + x^2 + 1*x + 1"))
factorize.S1P6.F3F3(toPoly.pm("x^6 + x^2 + 8*x + 1"))
sapply(1:7, function(x) factorize.S1P6.F3F3(toPoly.pm("x^6 + x^2 + x[1]*x + 1"))$isF)


############

### Factors: (x^3 + ... - 1)
# (x^3 + b2*x^2 + b1*x - 1) * (x^3 + c2*x^2 + c1*x - 1)
p = toPoly.pm("(x^3 + b2*x^2 + b1*x - 1) * (x^3 + c2*x^2 + c1*x - 1)");

b = c(3, -4); c = c(2, 5);
b1 = b[1]; b2 = b[2]; c1 = c[1]; c2 = c[2];
P = function(x) eval.pm(p, list(x=x, b1=b1, b2=b2, c1=c1, c2=c2));

#########
### Mod 7

### P(1) (mod 7)
((b2 + b1)*(c2 + c1) - P(1)) %% 7

### P(2) (mod 7)
(9*(3*b2 + b1)*(3*c2 + c1) - P(2)) %% 7
((3*b2 + b1)*(3*c2 + c1) + 3*P(2)) %% 7

### TODO


#################

### [old]
# - unsuccessful;

#########
### Mod 9
# - offers some interesting reductions: 2^3 + 1, 8^3 + 1;
# - BUT NO pathway to solve;

### P(2) (mod 9)
(4*(2*b2 + b1)*(2*c2 + c1) - P(2)) %% 9
((2*b2 + b1)*(2*c2 + c1) + 2*P(2)) %% 9

### P(8) (mod 9)
((b2 - b1)*(c2 - c1) - P(8)) %% 9

### P(1/2) (mod 9)
((-2*b2 - 4*b1)*(-2*c2 - 4*c1) - P(5)) %% 9
((b2 + 2*b1)*(c2 + 2*c1) + 2*P(5)) %% 9

### Aux: P(1) (mod 9)
((b2 + b1 + 2)*(c2 + c1 + 2) - P(1)) %% 9
### Aux: P(4) (mod 9)
(4*(b2 - 2*b1 - 1)*(c2 - 2*c1 - 1) - P(4)) %% 9
((b2 - 2*b1 - 1)*(c2 - 2*c1 - 1) + 2*P(4)) %% 9
### Aux: P(3) (mod 9)
((3*b1 + 1)*(3*c1 + 1) - P(3)) %% 9

# System =>
(b1*c1 + 4*b2*c2 + 2*(b2*c1 + b1*c2) + 2*P(2)) %% 9
(4*b1*c1 + b2*c2 + 2*(b2*c1 + b1*c2) + 2*P(5)) %% 9
(b1*c1 + b2*c2 - (b2*c1 + b1*c2) - P(8)) %% 9
(b1*c1 + b2*c2 + (b2*c1 + b1*c2) + 2*(b1+b2+c1+c2) + 4 - P(1)) %% 9
(4*b1*c1 + b2*c2 - 2*(b2*c1 + b1*c2) + (2*b1 - b2 + 2*c1 - c2) + 1 + 2*P(4)) %% 9
(3*(b1+c1) + 1 - P(3)) %% 9
# Linear System =>
(P(2) + P(5) - 2*P(8)) %% 9 # == 0!
# Sum: Eq 3, Eq 4 =>
(2*(b1*c1 + b2*c2) + 2*(b1+b2+c1+c2) + 4 - P(1) - P(8)) %% 9
((b1*c1 + b2*c2) + (b1+b2+c1+c2) + 2 + 4*P(1) + 4*P(8)) %% 9
# Eq  5 - Eq 6 =>
(4*b1*c1 + b2*c2 - 2*(b2*c1 + b1*c2) - (b1 + b2 + c1 + c2) + 2*P(4) + P(3)) %% 9
# Eq 1 + (Eq  5 - Eq 6) =>
(5*b1*c1 + 5*b2*c2 - (b1 + b2 + c1 + c2) + 2*P(2) + 2*P(4) + P(3)) %% 9
(2*b1*c1 + 2*b2*c2 + 2*(b1+b2+c1+c2) + 4 - P(1) - P(8)) %% 9
# =>
(3*(b1+b2+c1+c2) + 8 - 2*P(1) - 2*P(8) + 2*P(2) + 2*P(4) + P(3)) %% 9
(3*(b2+c2) + 7 - 2*P(1) - 2*P(8) + 2*P(2) + 2*P(4) + 2*P(3)) %% 9

#
(5*(b1*c1 + b2*c2) + 4*(b2*c1 + b1*c2) + 2*P(2) + 2*P(5)) %% 9
((b1*c1 + b2*c2) - (b2*c1 + b1*c2) + 4*P(2) + 4*P(5)) %% 9 # redundant


# Always True:
P2 = function(x, b1=-2) x^6 - b1*x + 1
(P2(2) + P2(5) - 2*P2(8)) %% 9

p2 = toPoly.pm("(x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1) * (x^2 - 5*x + 1)");
P2 = function(x) eval.pm(p2, list(x=x, b1=b1, b2=b2, c1=c1, c2=c2));
p3 = toPoly.pm("(x^8 + x^2 - 2*x + 1)");
P3 = function(x) eval.pm(p3, list(x=x));
# Always True:
(P2(2) + P2(5) - 2*P2(8)) %% 9
(P3(2) + P3(5) - 2*P3(8)) %% 9

