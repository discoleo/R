########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Modular Arithmetic
###
### draft v.0.1a


# - minimal Modular Arithmetic;
# - used to Factorize polynomials;


# - this file:
#   source("Polynomials.Helper.Mod.R");
# - is automatically loaded in:
#   Polynomials.Helper.Factorize.R;

##########################
##########################

### Inverse (mod p)
inv.mod = function(x, mod) {
	if(mod == 1 || mod == 0) stop("Invalid mod!");
	f = function(x) {
		id = which((x * seq(mod - 1) - 1) %% mod == 0);
		if(length(id) == 0) NA else id;
	}
	sapply(x, f);
}

### Solve P2:
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

