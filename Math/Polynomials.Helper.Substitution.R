########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Multi-Variable Polynomials
### Polynomial Substitution


### fast load:
# - is automatically loaded in Polynomials.Helper.R;
# source("Polynomials.Helper.Substitution.R")


#################################
#################################

### Replace/Substitute variables:

### Use Cases:
# - x => value;
# - m^pow => value; e.g. m^pow = 1,
#   where m = root of unity of order pow;
# - k^pow => "K", where k = K^(1/pow);
# - x => p2 (a polynomial);
### Shortcuts:
# - shift variable: x => (x + val);
# - scale variable: x => (val*x);
# - scale by root of unity: x => (m*x),
#   where m^pow = 1;


###########################

### Transforms

### Shift variables:
# x => (x + val)
shift.pm = function(p, val, xn="x", tol=1E-10) {
	len = length(val);
	if(len > 1) {
		if(length(xn) == 1) {
			warning("Same variable used!");
			return(shift.pm(p, round0(sum(val), tol=tol), xn=xn, tol=tol))
		}
		for(i in seq(len)) {
			p = shift.pm(p, val = val[i], xn=xn[i], tol=tol);
		}
		return(p);
	}
	# 1 Variable:
	if(val == 0) return(p);
	x.new = data.frame(x=1:0, coeff = c(1, val));
	names(x.new)[1] = xn;
	p = replace.pm(p, x.new, xn=xn, pow=1);
	if( ! is.null(tol)) p$coeff = round0(p$coeff, tol=tol);
	return(p);
}

# Note: scale (Generic) behaves badly;
rescale.pm = function(p, val, xn="x", mod=NULL, div=FALSE) {
	idx = match(xn, names(p));
	if(is.na(idx)) stop("Variable not found!");
	hasX = (p[, idx] != 0);
	if(div) {
		pow.max = maxPow.pm(p, xn);
		p$coeff = p$coeff * val^(pow.max - p[, idx]);
		if( ! is.null(mod)) stop("Mod NOT yet implemented!");
	} else {
		p$coeff[hasX] = p$coeff[hasX] * val^p[hasX, idx];
		if( ! is.null(mod)) {
			if(mod < 2) stop("Modulus must be >= 2!");
			p$coeff = p$coeff %% mod;
			p = reduce.pm(p);
		}
	}
	return(p);
}
# x => (x*m), where m^pow = 1;
rescale.pm.unity = function(p, pow, mn="m", xn="x", div=FALSE) {
	idx = match(xn, names(p));
	if(is.na(idx)) stop("Variable not found!");
	hasX = (p[, idx] != 0);
	x.pow = p[hasX, idx] %% pow;
	if(div) x.pow[x.pow != 0] = pow - x.pow[x.pow != 0];
	#
	idm = match(mn, names(p));
	if(is.na(idm)) {
		p[, mn] = 0;
		p[hasX, mn] = x.pow;
	} else {
		p[hasX, idm] = p[hasX, idm] + x.pow;
	}
	return(p);
}


# TODO:
# - explore:
#   methods::setGeneric("replace")
#   methods::setMethod("replace", signature = c(p1 = "pm", p2 = "numeric"), definition = replace.pm.numeric)
#   methods::setMethod("replace", signature = c(p1 = "pm", p2 = "character"), definition = replace.pm.character)
replace.withVal.pm = function(p, ...) stop("Defunct function: replace.withVal!");
replace.pm.numeric = function(p1, p2, xn, pow=1, simplify=TRUE, tol=1E-10) {
	p = p1; val = p2;
	if(missing(xn)) {
		xn = names(val);
		if(is.null(xn)) stop("Missing variable name!");
	}
	if(length(val) > 1) {
		len = length(val);
		if(length(pow) == 1) pow = rep(pow, len);
		if(length(xn) == 1) {xn = rep(xn, len); warning("Same variable used!");}
		for(i in seq(len)) {
			p = replace.pm.numeric(p, val[i], xn=xn[i], pow=pow[i], simplify=simplify, tol=tol);
		}
		return(p);
	}
	if(val == 0) {
		if(pow != 1) {
			warning("Only some terms will be replaced with 0!");
			isZero = p[,xn] >= pow; # every x^pow == 0;
			p = p[ ! isZero, , drop=FALSE];
		} else p = p[p[,xn] == 0, , drop=FALSE];
		return(p);
	}
	xpow = p[,xn];
	p[,xn] = if(pow == 1) 0 else xpow %% pow;
	hasX  = xpow >= pow; xpow = xpow[hasX];
	# Replace with value:
	x.pow = if(pow == 1) xpow else xpow %/% pow;
	# Optimized
	xpow.unq = sort(unique(x.pow));
	id = match(x.pow, xpow.unq);
	xval = val^xpow.unq;
	p[hasX, "coeff"] = p[hasX, "coeff"] * xval[id];
	p = aggregate0.pm(p);
	if(tol > 0) p$coeff = round0(p$coeff, tol=tol);
	p = reduce.pm(p);
	if(simplify) p = reduce.var.pm(p);
	return(p)
}
replace.pm = function(p1, p2, xn, pow=1, sequential=TRUE) {
	# replace x^pow by p2;
	if(is.numeric(p2) || is.complex(p2)) return(replace.pm.numeric(p1, p2=p2, xn=xn, pow=pow));
	if(is.character(p2)) return(replace.pm.character(p1, p2=p2, xn=xn, pow=pow, sequential=sequential));
	# Checks
	stop.f = function() stop("Missing variable name!");
	if(missing(xn)) {
		if(is.pm(p2)) stop.f();
		# Named List:
		if(is.list(p2)) {
			xn = names(p2);
			if(is.null(xn)) stop.f();
			len = length(p2);
			if(len > 1) {
				warning("More than 1 value: replacing sequentially!");
				if(length(pow) == 1) pow = rep(pow, len);
				for(id in seq(len)) {
					p1 = replace.pm(p1, p2[[id]], xn=xn[[id]], pow=pow[[id]]);
				}
				return(p1);
			}
			p2 = p2[[1]];
		} else stop("Replacement type: not supported!");
	}
	if(length(xn) > 1 || length(pow) > 1) stop("Only 1 value supported!")
	idx = match(xn, names(p1));
	if(any(is.na(idx))) {
		warning(paste0("Polynomial does NOT contain variable: ", xn));
		idx = idx[ ! is.na(idx)];
		if(length(idx) == 0) return(p1);
	}
	# xPow
	rpow = if(pow == 1) p1[,idx] else p1[,idx] %/% pow;
	p1[,idx] = if(pow == 1) 0 else p1[,idx] %% pow;
	max.pow = max(rpow);
	# Powers of p2:
	p2.pows = list(p2);
	if(max.pow > 1) {
		if(ncol(p2) > 1) {
			for(ipow in seq(2, max.pow)) {
				p2.pows[[ipow]] = mult.pm(p2.pows[[ipow - 1]], p2);
			}
		} else {
			p2.pows = sum(p2$coeff)^seq(max.pow);
		}
	}
	# Result:
	pR = data.frame();
	for(nr in seq(nrow(p1))) {
		pR = if(rpow[nr] == 0) add.pm(pR, p1[nr,])
			else add.pm(pR, mult.pm(p2.pows[[rpow[nr]]], p1[nr,]));
	}
	return(reduce.var.pm(pR));
}
# Replace variable with another Variable
replace.pm.character = function(p1, p2, xn, pow=1, sequential=TRUE, na.stop=TRUE) {
	len = length(p2);
	lenpow = length(pow); lenx = length(xn);
	if(len > 1) {
		if(lenx == 1 && lenpow > 1) {
			# x^seq() => replaced with c(...);
			# - may be an issue with sequential processing!
			# - unknown use cases;
			xn = rep(xn, lenpow); lenx = lenpow;
		}
		if(lenx == 1 && ! sequential) {
			stop("The same variable name!");
		} else if(lenx == 1 && sequential) {
			warning("The same variable name!");
			p2 = p2[p2 != xn];
			if(length(p2) >= 1) {
				p2 = p2[[1]]; len = 1;
			} else {
				warning("No replacement!");
				return(p1);
			}
		}
	} else if(lenx > 1) {
		p2 = rep(p2, lenx);
	}
	# Replacement:
	if(any(pow > 1))
		return(replaceByPow.pm.character(p1, p2, xn=xn, pow=pow, sequential=sequential));
	# same Power:
	return(replaceNames.pm(p1, p2, xn=xn, sequential=sequential));
}
### Replace simple monomial m with xn;
# m = character vector with the names of the variables;
# - similar to replace.pm.character, but pow=1;
# - TODO: handle overlap between m & xn;
replace.pm.m = function(p, m, xn, warn=TRUE, simplify=TRUE) {
	if(length(xn) != 1) stop("Replacement must be exactly 1 variable!");
	nms = names(p);
	idV = match(m, nms);
	if(any(is.na(idV))) {
		if(warn) warning("Missing variables ", m[is.na(idV)], "!");
		return(p);
	}
	#
	nr = nrow(p);
	minPow = sapply(seq(nr), function(nr) {
		min(p[nr, idV]);
	})
	hasM = (minPow != 0);
	if( ! any(hasM)) return(p);
	# Add new variable:
	if(is.na(match(xn, nms))) {
		p[ , xn] = 0;
	}
	p[ , xn] = p[ , xn] + minPow;
	# Update var-names:
	idV = match(m, names(p));
	for(id in idV) {
		p[ , id] = p[ , id] - minPow;
	}
	#
	if(simplify) {
		p = drop.pm(p);
		p = aggregate0.pm(p);
	}
	return(p);
}

### Rename vars in p1 with those in p2
# xn = old names;
replaceNames.pm = function(p1, p2, xn, sequential=FALSE, debug=TRUE) {
	len = length(xn);
	nms = dimnames.pm(p1);
	if(sequential) {
		for(id in seq(len)) {
			nms[nms == xn[id]] = p2[id];
		}
	} else {
		id = match(xn, nms);
		if(any(is.na(id))) {
			if(na.stop)
				stop(paste0("Variable names not found: ", xn[is.na(id)]));
			isVar = ! is.na(id);
			id = id[isVar];
			xn = xn[isVar];
			if(length(id) == 0) warning("No replacements!")
		}
		nms[id] = p2;
	}
	if(debug) print(paste0("New names: ", paste0(nms, collapse=", ")));
	# set new names:
	idCoeff = which(names(p1) == "coeff");
	names(p1)[ - idCoeff] = nms;
	### Duplicate names:
	isDuplic = duplicated(nms);
	if( ! any(isDuplic)) return(p1);
	# merge duplicates:
	duplNms = nms[isDuplic];
	ids = seq(ncol(p1)); ids = ids[ - idCoeff];
	id0 = ids; ids = ids[isDuplic];
	id0 = id0[match(duplNms, nms)];
	for(nc in seq_along(ids)) {
		p1[, id0[nc]] = p1[, id0[nc]] + p1[, ids[nc]];
	}
	p1 = p1[, ! isDuplic, drop=FALSE];
	p1 = aggregate0.pm(p1);
	p1 = reduce.pm(p1);
	return(toPoly.pm(p1));
}
# Replace: some powers != 1;
replaceByPow.pm.character = function(p1, p2, xn, pow=1, sequential=TRUE, reduce=TRUE, debug=TRUE) {
	if(sequential) {
		for(id in seq(length(xn))) {
			xold = xn[[id]];
			idx = match(xold, names(p1));
			if(is.na(idx)) {
				if(debug) print(paste0("Var name ", xold, " not found!"));
				next;
			}
			# Values:
			isPow1 = (pow[[id]] == 1);
			tmp = if(isPow1) p1[, idx] else p1[, idx] %/% pow[[id]];
			p1[, idx] = if(isPow1) 0 else p1[, idx] %% pow[[id]];
			# New Var
			xnew = p2[[id]];
			idnew = match(xnew, names(p1));
			if(is.na(idnew)) {
				p1[, xnew] = tmp;
			} else {
				p1[, idnew] = p1[, idnew] + tmp;
			}
		}
	} else {
		isDuplic = duplicated(xn);
		if(any(isDuplic)) stop(paste0("Variable names cannot be duplicated!", xn[isDuplic]));
		# Change vars:
		newDF = data.frame();
		for(id in seq(length(xn))) {
			# Initial Columns:
			xold = xn[[id]];
			idx = match(xold, names(p1));
			if(is.na(idx)) {
				if(debug) print(paste0("Var name ", xold, " not found!"));
				next;
			}
			# Values:
			isPow1 = (pow[[id]] == 1);
			tmp = if(isPow1) p1[, idx] else p1[, idx] %/% pow[[id]];
			p1[, idx] = if(isPow1) 0 else p1[, idx] %% pow[[id]];
			# New Var
			xnew = p2[[id]];
			idnew = match(xnew, names(newDF));
			if(is.na(idnew)) {
				newDF[, xnew] = tmp;
			} else {
				newDF[, idnew] = newDF[, idnew] + tmp;
			}
		}
		# merge newDF back into p1:
		xnew  = names(newDF); # should be unique
		idnew = match(xnew, names(p1));
		isName = ! is.na(idnew);
		# existing Vars:
		prevNames = idnew[isName];
		p1[, prevNames] = p1[, prevNames] + newDF[, prevNames];
		# new Vars:
		newNames = idnew[ ! isName];
		p1[, newNames] = newDF[, newNames];
	}
	if(reduce) p1 = reduce.pm(p1);
	return(p1);
}
# Replace: monomial pv with the variable pn;
replace.pm.character.pm = function(p1, pn, pv, ..., drop=TRUE) {
	if( ! is.pm(pv)) stop("The replaced monomial must be a polynomial!");
	if(nrow(pv) != 1) stop("Only monomials can be replaced!");
	# Names of pv:
	pv = drop.pm(pv);
	nms = names(pv);
	idc = match("coeff", nms);
	nms = nms[ - idc]; coeff = pv[, idc]; pv = pv[, - idc, drop=FALSE];
	if(length(nms) == 0) {
		warning("Nothing to replace!");
		return(p1); # TODO: replace numeric values;
	}
	### Names in p1
	id.nms = match(nms, names(p1));
	isNA = is.na(id.nms);
	if(any(isNA)) {
		warning(paste0("Not all variables present:\n", nms[isNA]));
		return(p1);
	}
	px = p1[, id.nms, drop=F];
	p1 = p1[, - id.nms, drop=F];
	### Powers:
	class(pv) = "data.frame"; # TODO: hack;
	isPow1 = (pv[1, ] == 1);
	hasPowGt1 = any( ! isPow1);
	if(hasPowGt1) {
		pP = sapply(seq(length(id.nms))[ ! isPow1], function(nc) {
			px[, nc] %/% pv[, nc];
		})
		min.pow = sapply(seq(nrow(px)), function(nr) {
			min(px[nr, isPow1], pP[nr, ]);
		})
	} else {
		min.pow = sapply(seq(nrow(px)), function(nr) {
			min(px[nr, ]);
		})
	}
	# New variable:
	id.new = match(pn, names(p1));
	if(is.na(id.new)) {
		p1[, pn] = min.pow;
	} else {
		p1[, id.new] = p1[, id.new] + min.pow;
	}
	# Adjust remaining Powers:
	for(nc in seq(length(id.nms))) {
		px[, nc] = px[, nc] - min.pow * pv[, nc];
	}
	# TODO: check name overlap between px & pn;
	p1 = cbind(p1, px);
	# Adjust coefficients:
	if(coeff != 1) {
		p1$coeff = p1$coeff / coeff^min.pow;
	}
	# Clean-up:
	if(drop) p1 = drop.pm(p1);
	return(p1)
}
# Replace: with fraction p2/p2fr
replace.fr.pm = function(p1, p2, p2fr, xn, pow=1) {
	# replace xn^pow by p2/p2fr;
	idx = match(xn, names(p1));
	if(is.na(idx)) stop(paste0("Polynomial does NOT contain variable: ", xn));
	# xPow
	rpow = if(pow == 1) p1[,idx] else p1[,idx] %/% pow;
	p1[,idx] = if(pow == 1) 0 else p1[,idx] %% pow;
	max.pow = max(rpow);
	# powers of p2 & p2fr:
	p2.pows = list(p2);
	p2fr.pows = list(p2fr); # powers of p2fr
	if(max.pow > 1) {
		for(ipow in seq(2, max.pow)) {
			print(paste0("Pow = ", ipow))
			p2.pows[[ipow]] = mult.pm(p2.pows[[ipow - 1]], p2);
			p2fr.pows[[ipow]] = mult.pm(p2fr.pows[[ipow - 1]], p2fr);
		}
		print("Finished Powers!");
	}
	p2m = c(tail(p2fr.pows, 1), p2.pows);
	print("Starting cross-multiplication:");
	# pre-align p1 & p2: verify if useful?
	pall = align.pm(p1, p2);
	p1 = pall[[1]];
	if(max.pow > 1) {
	for(ipow in seq(1, max.pow - 1)) {
		p2m[[ipow + 1]] = mult.pm(p2m[[ipow + 1]], p2fr.pows[[max.pow - ipow]]);
		p2m[[ipow + 1]] = align.pm(p1, p2m[[ipow + 1]], doReduce=FALSE)[[2]];
	}
	}
	print("Finished cross-multiplication!");
	pR = data.frame();
	for(nr in seq(nrow(p1))) {
		ipow = rpow[nr];
		tmp = mult.pm(p2m[[ipow + 1]], p1[nr,]);
		# TODO: separate sum
		pR = sum.pm(pR, tmp);
	};
	return(reduce.var.pm(pR));
}

