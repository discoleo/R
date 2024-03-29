########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Factorize: Derivations
###
### draft v.0.1j


### Factorize Multi-Variable Polynomials
### Various Derivations

### A. Univariate Polynomials


######################

### Helper Functions

source("Polynomials.Helper.R")
# - is automatically loaded in: Polynomials.Helper.R;
#   source("Polynomials.Helper.Factorize.R")


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

#############
### Mod 7 ###
#############

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
### Eqs:
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


### Factorize F3F3:
factorize.S1P6.F3F3 = function(p, debug=FALSE) {
	nms = names(p);
	idc = match("coeff", nms);
	if(is.na(idc)) stop("Not a polynomial!");
	xn = nms[ - idc];
	if(length(xn) != 1) stop("Not a univariate polynomial!");
	if(max(p[, xn]) != 6) stop("Not a polynomial of Order 6!");
	#
	negP = function(p) {
		isOdd = (p[, xn] %% 2 == 1);
		p$coeff[isOdd] = - p$coeff[isOdd];
		return(p);
	}
	mod = 7;
	r = factorize.S1P6.F3F3Mod7(p, debug=debug);
	if( ! r$isF) {
		p = negP(p);
		r = factorize.S1P6.F3F3Mod7(p, debug=debug);
		# NO Factors:
		if( ! r$isF) return(r);
		# Sign changes in the Factors:
		r$F[, c("b2", "c2")] = (mod - r$F[, c("b2", "c2")]);
		r$B0 = -1;
	} else {
		r$B0 = 1;
	}
	if( ! r$isF) return(r);
	### Mod 13;
	mod = 13;
	r2 = factorize.S1P6.F3F3Mod13(p, debug=debug);
	if( ! r2$isF) {
		if(r$B0 == -1) return(r2);
		# TODO: check also existence with B0 = -1?
		p = negP(p);
		r = factorize.S1P6.F3F3Mod7(p, debug=debug);
		if( ! r$isF) return(r2);
		print("TODO: Factor may still be possible!");
		return(r2);
	}
	if(r2$isF) {
		if(r$B0 == -1) {
			# revert Sign changes in the Factors:
			r2$F[, c("b2", "c2")] = (mod - r2$F[, c("b2", "c2")]);
		}
		r$F1 = r$F;
		r$F2 = r2$F; # TODO ???
		r$Mod = c(r$Mod, r2$Mod);
		r$F = solve.mod.S2P1.F3F3(r$F1, r$F2, mod=r$Mod);
		r$Mod = c(prod(r$Mod), r$Mod);
	}
	return(r);
}
factorize.S1P6.F3F3Mod7 = function(p, debug=FALSE) {
	mod = 7;
	p1x = eval.pm(p, 1); p3x = eval.pm(p, 3);
	p5x = eval.pm(p, 5); p6x = eval.pm(p, 6);
	p1 = p1x %% mod; p3 = p3x %% mod;
	p5 = p5x %% mod; p6 = p6x %% mod;
	b1x = (2 - 4*p1 + 2*p3 + 2*p5 - p6) %% mod;
	b0x = (4*p3 + 4*p5 - 2*p6) %% mod;
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
		printVars.V4();
	}
	### Solve Sub-System
	### & Test Aux Eqs:
	p2 = eval.pm(p, 2) %% mod;
	p4 = eval.pm(p, 4) %% mod;
	sol = solve.F3.Coeffs(b1c1, b2c2, b12s, c12s, mod=mod,
		test=list(x=c(2, 4), p=c(p2, p4)), debug=debug);
	return(sol);
}
factorize.S1P6.F3F3Mod13 = function(p, debug=FALSE) {
	mod = 13;
	p1x = eval.pm(p, 1); p12x = eval.pm(p, -1);
	p4x = eval.pm(p, 4); p10x = eval.pm(p, 10);
	p1 = p1x %% mod; p12 = p12x %% mod;
	p4 = p4x %% mod; p10 = p10x %% mod;
	# b1x = - (b12s + c12s);
	b1x = (2 - 7*p1 - 2*p12 + 4*p4 + 4*p10) %% mod;
	b0x = ((9*p12 + 8*p4 + 8*p10)) %% mod;
	#
	sol = solve.ModP2(c(b0x, b1x, 1), mod=mod);
	if( ! sol$hasSol) return(list(isF = FALSE, Mod = mod));
	b12s = sol$Sol;
	c12s = (- b1x - b12s) %% mod;
	b1c1 = (- 4*p12 + 3*p4 + p10) %% mod;
	b2c2 = (9*p12 + p4 + 3*p10) %% mod;
	len  = length(b12s);
	b1c1 = rep(b1c1, len); b2c2 = rep(b2c2, len);
	if(debug) {
		printVars.V4();
	}
	### Solve Sub-System
	### & Test Aux Eqs:
	vals = c(2,3);
	px = sapply(vals, function(x) eval.pm(p, x)) %% mod;
	sol = solve.F3.Coeffs(b1c1, b2c2, b12s, c12s, mod=mod,
		test=list(x=vals, p=px), debug=debug);
	return(sol);
}
solve.F3.Coeffs = function(b1c1, b2c2, b12s, c12s, mod, testVals=NULL, debug=FALSE) {
	# (b1 - b12s)*(c1 - c12s) - b2c2 = 0
	# c12s*b1^2 - (b1c1 - b2c2 + b12s*c12s)*b1 + b1c1*b12s = 0
	b1y = (- (b1c1 - b2c2 + b12s*c12s)) %% mod;
	b0y = (b1c1*b12s) %% mod;
	len = length(b0y);
	sol = lapply(seq(len), function(id) solve.ModP2(c(b0y[id], b1y[id], c12s[id]), mod=mod));
	# Filter Solutions:
	len = length(sol);
	b1  = unlist(lapply(sol, function(sol) if(sol$hasSol) sol$Sol else NULL));
	if(length(b1) == 0) return(list(isF = FALSE, Mod = mod));
	#
	nSols = unlist(lapply(seq(len),
		function(id) if(sol[[id]]$hasSol) rep(id, length(sol[[id]]$Sol)) else c()));
	# print(nSols); print(b1);
	b12s = b12s[nSols]; c12s = c12s[nSols];
	b1c1 = b1c1[nSols]; b2c2 = b2c2[nSols];
	b2 = (b12s - b1) %% mod;
	# c1 & c2:
	len = length(b1);
	c1 = rep(NA, len); c2 = rep(NA, len);
	# TODO: fails when both b1 = 0 & b2 = 0;
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
	### Aux Eqs:
	if( ! is.null(testVals) && nrow(sol) > 0) {
		if(debug) {
			print(paste0("Starting Aux eqs: sol = ", nrow(sol)));
			print(sol);
		}
		x = testVals$x; p = testVals$p;
		b1 = sol[,1]; b2 = sol[,2]; c1 = sol[,3]; c2 = sol[,4];
		isSol = TRUE;
		for(id in seq_along(x)) {
			xi = x[id];
			x3 = (xi^3 + 1) %% mod;
			x2 = xi^2;
			err = ((b2*x2 + b1*xi + x3)*(c2*x2 + c1*xi + x3) - p[id]) %% mod;
			isSol = isSol & (err == 0);
		}
		sol  = sol[ isSol, , drop=FALSE];
	}
	if(nrow(sol) == 0) return(list(isF=FALSE, Mod=mod));
	return(list(isF=TRUE, F=sol, Mod=mod));
}
### Other:
solve.mod.S2P1.F3F3 = function(x, x2, mod, simplify=TRUE) {
	if(length(mod) != 2) stop("Both primes are needed!");
	if(ncol(x) != ncol(x2)) stop("Number of variables must match!");
	if(nrow(x) > 1 && simplify) x = unique.matrix.perm(x, c(2,2));
	#
	gr = expand.grid(lapply(c(nrow(x), nrow(x2), ncol(x2)), seq));
	sol = sapply(seq(nrow(gr)), function(id) {
		nc = gr[id, 3];
		vals = c(x[gr[id, 1], nc], x2[gr[id, 2], nc]);
		solve.ModP1Base(vals, mod=mod);
	});
	sol = matrix(unlist(sol), ncol=ncol(x2));
	# Filter 0-solution:
	pr = prod(mod);
	isZero = apply(sol, 1, function(x) all((x %% pr) == 0));
	sol = sol[ ! isZero, ];
	return(sol);
}
seq.tokens.perm = function(tokens) {
	if(length(tokens) == 1) return(seq(tokens));
	tkStart = cumsum(c(1, head(tokens, -1)));
	tkEnd   = cumsum(tokens);
	idTk = seq(length(tokens));
	idTk = c(idTk[-1], idTk[1]);
	tkStart = tkStart[idTk]; tkEnd = tkEnd[idTk];
	tk = unlist(lapply(seq(along=tkStart), function(id) seq(tkStart[id], tkEnd[id])));
	return(tk);
}
unique.matrix.perm = function(x, tokens=NULL) {
	nr = nrow(x);
	if(nr == 1) return(x);
	if(nr == 2) {
		if(all(x[1,] == x[2,])) return(x[1, , drop=FALSE]);
		if(is.null(tokens)) return(x);
		# Permutation:
		# TODO: more than 2 tokens to permute;
		id  = seq.tokens.perm(tokens);
		tmp = rbind(x, x[ , id]);
		isDuplicated = duplicated(tmp);
		if(any(isDuplicated)) return(x[1, , drop=FALSE]);
		return(x);
	}
	# Case: nr > 2
	x   = unique(x);
	nr0 = nrow(x);
	if(nr0 == 1) return(x);
	# Check permutations:
	idc = seq.tokens.perm(tokens);
	isSol = rep(TRUE, nr0);
	for(nr in seq(nr0, 2)) {
		tmp  = x[nr, idc];
		isEq = sapply(seq(nr - 1), function(nr2) all(x[nr2,] == tmp));
		if(any(isEq)) isSol[nr] = FALSE;
	}
	x = x[isSol, , drop=FALSE];
	return(x)
}
### Debug:
printVars.V4 = function() {
	# debug various functions;
	l = parent.frame();
	with(l, {
	cat(c("\nEq 1: ", c(b1x, b0x)));
	cat(c("\nb1 + b2 = ", b12s));
	cat(c("\nc1 + c2 = ", c12s));
	cat(c("\nb1 * c1 = ", b1c1));
	cat(c("\nb2 * c2 = ", b2c2, "\n"));
	});
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

### False-Positives: with (mod 7)
# What is the meaning?
# - much less frequent with both (mod 7) & (mod 13);
factorize.S1P6.F3F3(toPoly.pm("x^6 + x^2 + 1*x + 1"))
factorize.S1P6.F3F3(toPoly.pm("x^6 + x^2 + 8*x + 1"))
sapply(1:22, function(x) factorize.S1P6.F3F3(toPoly.pm("x^6 + x^2 + x[1]*x + 1"))$isF)


############

### Factors: (x^3 + ... - 1)
# (x^3 + b2*x^2 + b1*x - 1) * (x^3 + c2*x^2 + c1*x - 1)
p = toPoly.pm("(x^3 + b2*x^2 + b1*x - 1) * (x^3 + c2*x^2 + c1*x - 1)");

b = c(3, -4); c = c(2, 1);
b1 = b[1]; b2 = b[2]; c1 = c[1]; c2 = c[2];
p2 = replace.pm(p, list(b1=b1, b2=b2, c1=c1, c2=c2))
P = function(x) eval.pm(p2, x);
factorize.S1P6.F3F3(p2)

###
b = c(3, -1); c = c(3, 1);
b1 = b[1]; b2 = b[2]; c1 = c[1]; c2 = c[2];
p2 = replace.pm(p, list(b1=b1, b2=b2, c1=c1, c2=c2))
P = function(x) eval.pm(p2, x);
factorize.S1P6.F3F3(p2)

# NO False-Positives:
sapply(1:7, function(b) factorize.S1P6.F3F3(toPoly.pm("x^6 + b[1]*x + 1"))$isF)
# False-Positives: ~ 1/4 of tests with (mod 7);
sapply(1:100, function(b) factorize.S1P6.F3F3(toPoly.pm("x^6 + x^2 + b[1]*x + 1"))$isF)


# Note:
# - can be solved using the same code as for:
#   (x^3 + ... + 1)*(x^3 + ... + 1);
# - Transform: P(-x);

#########
### Mod 7

### P(1) (mod 7)
((b2 + b1)*(c2 + c1) - P(1)) %% 7
(b1*c1 + b2*c2 + (b1*c2 + b2*c1) - P(1)) %% 7

### P(2) (mod 7)
(4*(2*b2 + b1)*(2*c2 + c1) - P(2)) %% 7
((2*b2 + b1)*(2*c2 + c1) - 2*P(2)) %% 7
(b1*c1 + 4*b2*c2 + 2*(b1*c2 + b2*c1) - 2*P(2)) %% 7

### P(4) (mod 7)
(2*(4*b2 + b1)*(4*c2 + c1) - P(4)) %% 7
((4*b2 + b1)*(4*c2 + c1) - 4*P(4)) %% 7
(b1*c1 + 2*b2*c2 + 4*(b1*c2 + b2*c1) - 4*P(4)) %% 7

# System =>
(b1*c1 - 2*b2*c2 - 2*P(1) + 2*P(2)) %% 7
(3*b1*c1 + 2*b2*c2 - 4*P(1) + 4*P(4)) %% 7

### Eqs:
(b1*c1 + 2*P(1) + 4*P(2) + P(4)) %% 7
(b2*c2 + 2*P(1) - 6*P(2) + 4*P(4)) %% 7
((b1*c2 + b2*c1) - 5*P(1) + 2*P(2) - 5*P(4)) %% 7

# - solution using P(-x);


#################
#################

##############
### Mod 13 ###
##############

### Factors: (x^3 + ... + 1)
# (x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1)
p = toPoly.pm("(x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1)");

b = c(3, -4); c = c(2, 1);
b1 = b[1]; b2 = b[2]; c1 = c[1]; c2 = c[2];
p2 = replace.pm(p, list(b1=b1, b2=b2, c1=c1, c2=c2))
P = function(x) eval.pm(p2, x);

### P(12) or P(-1) (mod 13)
((b2 - b1)*(c2 - c1) - P(-1)) %% 13
(b1*c1 + b2*c2 - (b1*c2 + b2*c1) - P(-1)) %% 13

### P(4) (mod 13)
(3*(4*b2 + b1)*(4*c2 + c1) - P(4)) %% 13
((4*b2 + b1)*(4*c2 + c1) + 4*P(4)) %% 13
(b1*c1 + 3*b2*c2 + 4*(b1*c2 + b2*c1) + 4*P(4)) %% 13

### P(10) (mod 13)
(9*(10*b2 + b1)*(10*c2 + c1) - P(10)) %% 13
((10*b2 + b1)*(10*c2 + c1) - 3*P(10)) %% 13
(b1*c1 + 9*b2*c2 - 3*(b1*c2 + b2*c1) - 3*P(10)) %% 13

# System =>
(5*b1*c1 + 7*b2*c2 - 4*P(-1) + 4*P(4)) %% 13
(b1*c1 - 3*b2*c2 + 5*P(-1) - 5*P(10)) %% 13

### Eqs:
(b1*c1 + 4*P(-1) - 3*P(4) - P(10)) %% 13
(b2*c2 - 9*P(-1) - P(4) - 3*P(10)) %% 13
((b1*c2 + b2*c1) - 4*P(-1) - 4*P(4) - 4*P(10)) %% 13

### Aux: P(1) (mod 13)
((b2 + b1 + 2)*(c2 + c1 + 2) - P(1)) %% 13
((b1+b2)*(c1+c2) + 2*(b1+b2+c1+c2) + 4 - P(1)) %% 13
# =>
((b1+b2)*(c1+c2) - 9*P(-1) - 8*P(4) - 8*P(10)) %% 13
((b1+b2+c1+c2) + 2 - 7*P(1) - 2*P(-1) + 4*P(4) + 4*P(10)) %% 13
# =>
x = (b1 + b2)
(x^2 + (2 - 7*P(1) - 2*P(-1) + 4*P(4) + 4*P(10))*x + (9*P(-1) + 8*P(4) + 8*P(10))) %% 13


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

