#####################
##
## Leonard Mada
## (the one and only)
##
## DE: PDE
## Basic Types
##
## v.0.1a


####################

### Helper Functions


col = c("#FAB664A0", "#FAF664A0", "#B2FA64A0")

### Differentiate
D.expr = function(x, by) {
	y = substitute(expression(x));
	y = y[[2]];
	dx = lapply(by, function(nm) {
		D(y, nm);
	})
	return(dx[[1]])
}

###################

### Diffusion:
# Pt = D/2 * Pxx

# D2 = 2*D;
Pt  = D.expr(exp(-x^2 / (D2*t)) / sqrt(D2*t), "t")
Px  = D.expr(exp(-x^2 / (D2*t)) / sqrt(D2*t), "x")
Pxx = D(Px, "x")

### Test:
# Pt = D/2 * Pxx
x = sqrt(3); D = sqrt(5)/2; t = sqrt(pi); D2 = 2*D; # just a Test;
params = list(x=x, D2 = 2*D, t=t);
eval(Pt, params)
eval(Pxx, params) * D / 2;


# Symbolic Derivation:

# Pt:
exp(-x^2/(D2 * t)) * (x^2 * D2/(D2 * t)^2)/sqrt(D2 * t) +
	- exp(-x^2/(D2 * t)) * (0.5 * (D2 * (D2 * t)^-0.5))/sqrt(D2 * t)^2;
D2/2 * exp(-x^2/(D2 * t)) * (2*x^2 - D2*t)/(D2 * t)^2 / sqrt(D2 * t);
D/2 * exp(-x^2/(2*D * t)) * (x^2 - D*t)/(D * t)^2 / sqrt(2*D * t);

# Pxx:
eval(Pt, params) * 2 / D;
-(exp(-x^2/(D2 * t)) * (2/(D2 * t)) - exp(-x^2/(D2 * t)) *
	(2 * x/(D2 * t)) * (2 * x/(D2 * t)))/sqrt(D2 * t);
2*exp(-x^2/(D2 * t)) * (2*x^2 - D2*t)/(D2 * t)^2 / sqrt(D2 * t);
exp(-x^2/(2*D * t)) * (x^2 - D*t)/(D * t)^2 / sqrt(2*D * t);

