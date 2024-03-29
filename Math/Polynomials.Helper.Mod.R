########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Modular Arithmetic
###
### draft v.0.1g-root2-enh2


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


filter.mod = function(x, r, mod, exclude=FALSE) {
	len = length(mod);
	if(len == 1) {
		r0 = x %% mod;
		if(length(r) == 1) {
			isMod = (r0 == r);
		} else {
			isMod = (r0 %in% r);
		}
		if(exclude) isMod = ! isMod;
		return(x[isMod]);
	}
	# multiple Congruences
	# TODO
}

# Multiple:
# - each equation x^pow = valid y,
#   has multiple solutions;
# Strict All: x^pow = r, solvable for any r;
# Mixed: type All in one factor and type Multiple in another factor;
primes.mod = function(pow, type="Multiple", to=1000) {
	type = pmatch(type, c("Any", "Multiple", "Strict All", "Mixed"));
	if(is.na(type)) stop("Type is not supported!");
	p = primes(to); # library pracma;
	#
	if(type == 1) return(p);
	#
	warn.f = function(p) {
		warning("Mixed: not valid for power = strict prime!");
		# TODO: return "Multiple" or return c();
		p = numeric(0);
		attr(p, "pow") = pow;
		return(p);
	}
	if(pow == 3) {
		isMod = (p %% 6) != 5;
		if(type == 4) return(warn.f(p[isMod]));
	} else if(pow %% 2 == 0) {
		pFact = factors(pow);
		p2    = (pFact == 2);
		pFact = pFact[ ! p2];
		p2 = sum(p2);
		if(length(pFact) == 0) {
			isMod = (p %% (2^p2) == 1); # Multiple
			# Mixed:
			if(type == 4) { isMod = (p > 2) & ( ! isMod); }
			else if(type == 3) { p = 2; isMod = FALSE; }
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
			if(type == 4) return(warn.f(p[isMod]));
			if(type > 2) isMod = ! isMod;
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
			# Mixed: between All & Multiple values solvable;
			isMod = (p %% pFact[1]) == 1;
			isModM = isMod;
			for(id in seq(2, len)) {
				tmp = (p %% pFact[id] == 1);
				isMod = isMod | tmp;
				isModM = isModM & tmp;
			}
			isMod = isMod & ( ! isModM);
		}
		p = p[isMod];
		attr(p, "pow") = pow;
		return(p);
	}
	if(type == 2 || type == 4) { p = p[isMod]; }
	else if(type == 3) { p = p[ ! isMod]; }
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

### Solve S2P1:
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
			x = unique(x); # relevant for mod = c(2);
		} else {
			err = sapply(seq(mod), x.f);
			x   = which(err == 0);
			if(length(x) == 0) return(list(hasSol = FALSE, Mod = mod));
			x = (x + sh) %% mod;
		}
	}
	return(list(hasSol = TRUE, Sol = x, Mod = mod));
}

### Roots

### Order 2
root2.mod = function(x, mod) {
	# assumes: mod = prime;
	if(mod == 2) return(x);
	if(mod == 3) {
		r = (x %% mod);
		# TODO: vector;
		if(r == 2) return(NA);
		if(r == 1) return(c(1,2));
		return(0);
	}
	# Type of prime:
	rr = mod %% 4;
	if(rr == 3) {
		r = root2ModBase(x, mod=mod);
		r2 = (r*r) %% mod;
		if(r2[1] == x) return(r);
		return(NA);
	}
	### 1 (mod 4)
	rr = mod %% 8;
	if(rr == 5) {
		k = (mod + 3)/8;
		r = pow.mod(x, k, mod=mod);
		r2 = (r*r) %% mod;
		if(r2 == x) {
			# OK
		} else if(r2 + x == mod) {
			# TODO: use sqrt(-1);
			rn = pow.mod(2, 2*k-1, mod=mod);
			r  = (r*rn) %% mod;
		} else {
			return(NA); # NOT a quadratic residue!
		}
		r  = c(r, mod-r);
		return(r);
	}
	### 1 (mod 8)
	stop("Not yet implemented!");
	# TODO
}
root2ModBase = function(x, mod) {
	# assumes x is a quadratic residue!
	k = (mod + 1) / 4;
	if(length(x) == 1) {
		r = pow.mod(x, k, mod=mod);
		r = c(r, mod - r);
		if(r[2] < r[1]) r = c(r[2], r[1]);
	} else {
		r = lapply(x, function(x) {
			r = pow.mod(x, k, mod=mod);
			r = c(r, mod - r);
			if(r[2] < r[1]) r = c(r[2], r[1]);
			return(r)
		});
		r = do.call(rbind, r);
	}
	return(r);
}
### Order 3
root3.mod = function(x, mod) {
	# assumes: mod = prime;
	if(mod == 2 || mod == 3) return(x);
	# Type of prime:
	r6 = mod %% 6;
	if(r6 != 5) {
		stop("Not yet implemented!");
		# TODO
	}
	#
	return(root3ModBase(x, mod=mod));
}
root3ModBase = function(x, mod) {
	k   = (mod + 1) / 3;
	pow = mod - k;
	if(length(x) == 1) return(pow.mod(x, pow, mod=mod));
	r = sapply(x, pow.mod, pow, mod=mod);
	return(r);
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

####################
####################

