########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Modular Arithmetic
###
### draft v.0.1e


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
	if(x == 1) return(1);
	f = function(x) {
		id = which((x * seq(mod - 1) - 1) %% mod == 0);
		if(length(id) == 0) NA else id;
	}
	sapply(x, f);
}
inv2.mod = function(mod, pow=1) {
	if(mod %% 2 == 0) return(NA);
	xinv = (mod + 1) %/% 2;
	if(pow > 1) xinv = pow.mod(xinv, n=pow, mod=mod);
	return(xinv);
}
inv3.mod = function(mod, pow=1) {
	r = (mod %% 3);
	if(r == 0) return(NA);
	if(r == 2) {
		xinv = (mod + 1) %/% 3;
	} else xinv = mod - ((mod - 1) %/% 3);
	if(pow > 1) xinv = pow.mod(xinv, n=pow, mod=mod);
	return(xinv);
}

### Solve P2:
solve.ModP1 = function(r, mod) {
	# TODO: check gcd == 1;
	return(solve.ModP1Base(r, mod=mod));
}
solve.ModP1Base = function(r, mod) {
	len = length(r);
	if(len == 2) {
		rsol = ((r[1] - r[2]) * inv.mod(mod[2], mod[1])) %% mod[1];
		prM  = mod[1] * mod[2];
		rsol = (mod[2]*rsol + r[2]) %% prM;
		return(list(r=rsol, Mod=prM));
	} else if(len == 1) {
		r = r %% mod;
		return(list(r=r, Mod=mod));
	} else if(len > 2) {
		prM = mod[1]; rsol = r[1];
		for(id in seq(2, len)) {
			rsol = ((rsol - r[id]) * inv.mod(mod[id], prM)) %% prM;
			prM  = prM * mod[id];
			rsol = (mod[id]*rsol + r[id]) %% prM;
		}
		return(list(r=rsol, Mod=prM));
	}
	return(list(r=0, Mod=NA));
}

### Solve P2:
solve.ModP2 = function(b, mod) {
	b0 = b[1]; b1 = b[2]; b2 = b[3];
	# x.f = function(x) (b2*x^2 + b1*x + b0) %% mod;
	x.f = function(x) (x^2 + b0) %% mod;
	if(b0 == 0) {
		if(b2 == 1) {
			x = c(0, (mod - b1));
		} else if(b2 == 0) {
			x = 0;
		} else {
			b2inv = inv.mod(b2, mod=mod);
			x = c(0, (- b1*b2inv) %% mod);
		}
	} else if(b2 == 0) {
		x = - b0 * inv.mod(b1, mod=mod);
		x = x %% mod;
	} else {
		# Optimized version:
		if(b1 == 0) {
			sh = 0;
		} else {
			if(b2 != 1) {
				b2inv = inv.mod(b2, mod=mod);
				b1 = (b1 * b2inv) %% mod;
				b0 = (b0 * b2inv) %% mod;
			}
			inv2 = inv2.mod(mod, pow=1);
			sh = (- b1 * inv2) %% mod;
			b0 = (b0 - sh*sh) %% mod;
		}
		if((b0 + 1) %% mod == 0) {
			x = (1 + sh) %% mod;
		} else {
			err = sapply(seq(mod), x.f);
			x   = which(err == 0);
			if(length(x) == 0) return(list(hasSol = FALSE, Mod = mod));
			x = (x + sh) %% mod;
		}
	}
	return(list(hasSol = TRUE, Sol = x, Mod = mod));
}


####################

whichZero.mod = function(n, pow, b0=1, mod=n) {
	if(n > 1000) {
		x = seq(n);
		x = sapply(x, function(x) pow.mod(x, n=pow, mod=mod));
		x = (x + b0) %% mod;
	} else {
		x = (seq(n)^pow + b0) %% mod;
	}
	which(x == 0);
}
whichZero.bigz = function(n, pow, b0=1, mod=n) {
	# very slow for large numbers!
	x = gmp::as.bigz(seq(n), mod);
	x = lapply(x, function(x) pow.bigz(x, n=pow, mod=mod));
	x = do.call(c, x);
	x = x + b0;
	which(x == 0);
}
pow.mod = function(x, n, mod) {
	if(n == 1) return(x %% mod);
	if(is.double(n) && (n == round(n))) n = as.integer(n);
	if( ! is.integer(n)) stop("n must be integer!");
	# Multiply
	p.r   = 1;
	p.pow = x %% mod;
	while (n > 0) {
		if (n %% 2 == 1) {
			p.r = (p.r * p.pow) %% mod;
		}
		if(n == 1) break;
		p.pow = (p.pow * p.pow) %% mod;
		n = n %/% 2;
    }
	return(p.r);
}
pow.bigz = function(x, n, mod=NULL) {
	if( ! gmp::is.bigz(x)) {
		if(is.null(mod)) {
			x = gmp::as.bigz(x);
		} else {
			x = gmp::as.bigz(x, mod);
		}
	}
	gmp::pow.bigz(x, n);
}

