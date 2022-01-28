########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Modular Arithmetic
###
### draft v.0.1f-fix


# - minimal Modular Arithmetic;
# - used to Factorize polynomials;


# - this file:
#   source("Polynomials.Helper.Mod.R");
# - is automatically loaded in:
#   Polynomials.Helper.Factorize.R;
# - top file:
#   source("Polynomials.Helper.R");

##########################
##########################


filter.mod = function(x, r, mod) {
	len = length(mod);
	if(len == 1) {
		r0 = x %% mod;
		if(length(r) == 1) {
			isMod = (r0 == r);
		} else {
			isMod = (r0 %in% r);
		}
		return(x[isMod]);
	}
	# multiple Congruences
	# TODO
}
primes.mod = function(pow, type="Multiple", to=1000) {
	# Multiple = multiple solutions for each x^pow = valid solution;
	# All values: x^pow = r, solvable for any r;
	type = pmatch(type, c("Any", "Multiple", "Strict All", "Most/All values"));
	if(is.na(type)) stop("Type is not supported!");
	p = primes(to); # library pracma;
	#
	if(type == 1) return(p);
	if(pow == 3) {
		isMod = (p %% 6) != 5;
	} else if(pow %% 2 == 0) {
		pFact = factors(pow);
		p2    = (pFact == 2);
		pFact = pFact[ ! p2];
		p2 = sum(p2);
		if(length(pFact) == 0) {
			isMod = (p %% (2^p2) != 1);
		} else {
			# TODO
			stop("Not yet implemented!");
		}
	} else {
		pFact = factors(pow);
		pFact = unique(pFact); # TODO
		len   = length(pFact);
		if(len == 1) {
			isMod = (p %% (2*pow)) == 1;
		} else if(type == 3) {
			# Strict: All values solvable;
			isMod = (p %% pFact[1]) != 1;
			for(id in seq(2, len)) {
				isMod = isMod & (p %% pFact[id] != 1);
			}
		} else if(type == 2) {
			# Strict: few values solvable;
			isMod = (p %% pFact[1]) == 1;
			for(id in seq(2, len)) {
				isMod = isMod & (p %% pFact[id] == 1);
			}
		} else {
			# Few/In-between values solvable;
			isMod = (p %% pFact[1]) == 1;
			for(id in seq(2, len)) {
				isMod = isMod | (p %% pFact[id] == 1);
			}
		}
		p = p[isMod];
		attr(p, "pow") = pow;
		return(p);
	}
	if(type == 2) { p = p[isMod]; }
	else if(type == 3 || type == 4) { p = p[ ! isMod]; }
	else { return(numeric()); }
	attr(p, "pow") = pow;
	return(p);
}
countSol.mod = function(p, pow, type="Values") {
	type = pmatch(type, c("Count", "Diff", "Values"));
	if(is.na(type)) stop("Type NOT supported!");
	if(missing(pow)) {
		pow = attr(p, "pow");
		if(is.null(pow)) stop("pow is missing!")
	}
	#
	if(type == 1) {
		sapply(p, function(p) length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )));
	} else if(type == 2) {
		sapply(p, function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )));
	} else {
		rbind(p,
			sapply(p, function(p) length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
		)
	}
}

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
		if(b2 != 1) {
			b2inv = inv.mod(b2, mod=mod);
			b1 = (b1 * b2inv) %% mod;
			b0 = (b0 * b2inv) %% mod;
		}
		# Shift:
		if(b1 == 0) {
			sh = 0;
		} else {
			inv2 = inv2.mod(mod, pow=1);
			sh = (- b1 * inv2) %% mod;
			b0 = (b0 - sh*sh) %% mod;
		}
		if((b0 + 1) %% mod == 0) {
			x = c(sh + 1, sh - 1) %% mod;
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

