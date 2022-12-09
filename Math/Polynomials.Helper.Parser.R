########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Parse Polynomials
###
### draft v.0.2d

### Parser for Multi-variable Polynomials


##############

### fast load:
# Note:
# - is automatically loaded by: Polynomials.Helper.R;
# source("Polynomials.Helper.Parser.R")

###############

###############
### History ###
###############


### draft v.0.2c:
# - better handling of environments;
### draft v.0.2a - v.0.2a-ref:
# - moved from Polynomials.Helper.R;
# - plan for refactoring;
### draft v.0.1 [branch]
# - developed during:
#   Polynomials.Helper.R;


##########################
##########################

##############
### Parser ###
##############

### Parse expressions / polynomials
as.pm = function (p, ...) UseMethod("as.pm");
as.pm.default = function(p, ...) {
	# if(inherits(p, "polynomial")) {
		# return(as.pm.polynomial(p, ...));
	# }
	return(toPoly.pm(p, ...));
}

# Class polynomial => "pm"
# - from package polynom;
as.pm.polynomial = function(p, xn="x", sort=NULL, tol=1E-8) {
	if(length(xn) != 1) stop("Only univariate polynomials are supported!");
	p = unclass(p);
	if(tol != 0) p = round0(p, tol=tol);
	len = length(p);
	pR  = data.frame(x = seq(0, len-1), coeff = p);
	names(pR)[1] = xn;
	if( ! is.null(sort)) pR = sort.pm(pR, xn);
	class(pR) = c("pm", class(pR));
	return(pR);
}

# Read from Clipboard
polyClip = function(env=parent.frame(), reduce=FALSE, verbose=FALSE) {
	p = paste0(readClipboard(), collapse="\n");
	toPoly.pm(p, env=env, reduce=reduce, verbose=verbose);
}

### Parser:
toPoly.pm = function(e, env=NULL, reduce=FALSE, verbose=TRUE) {
	if(is.null(env)) env = parent.frame();
	if(is.character(e)) {
		if(length(e) > 1) {
			# List of polynomials
			pl = lapply(e, function(e) toPoly.pm(e, env=env));
			return(pl);
		}
		e = parse(text=e);
	} else if(inherits(e, "formula")) {
		e = e[[2]];
	} else if(is.numeric(e) || is.complex(e)) {
		p = if(length(e) == 1) data.frame(coeff=e)
			else stop("Not yet implemented!");
		class(p) = c("pm", class(p));
		return(p);
	} else if(inherits(e, "data.frame")) {
		p = e;
		if(reduce) p = reduce0.pm(p);
		if(inherits(p, "pm")) return(p);
		class(p) = c("pm", class(p));
		return(p);
	}
	if(is.expression(e)) {
		e = e[[1]];
	} else if( ! (is.language(e) || is.numeric(e) || is.complex(e)) ) {
		stop("Input must be an expression!");
	}
	p = data.frame();
	while(TRUE) {
		isSymbol = is.symbol(e);
		if(isSymbol || is.symbol(e[[1]])) {
			op = if(isSymbol) e else e[[1]];
			if(op == "+") {
				m = toMonom.pm(e[[3]], env=env, verbose=verbose);
				p = if(nrow(p) == 0) m else sum.pm(p, m);
				e = e[[2]];
			} else if(op == "-") {
				if(length(e) > 2) {
					m = toMonom.pm(e[[3]], xsign=-1, env=env, verbose=verbose);
					p = if(nrow(p) == 0) m else sum.pm(p, m);
					e = e[[2]];
				} else {
					m = toMonom.pm(e[[2]], xsign=-1, env=env, verbose=verbose);
					p = if(nrow(p) == 0) m else sum.pm(p, m);
					break;
				}
			} else {
				m = toMonom.pm(e, env=env, verbose=verbose);
				p = if(nrow(p) == 0) m else sum.pm(p, m);
				break;
			}
		} else if(is.numeric(e) || is.complex(e)) {
			m = toMonom.pm(e, env=env, verbose=verbose);
			p = if(nrow(p) == 0) m else sum.pm(p, m);
			break;
		} else break;
	}
	if(reduce) p = reduce0.pm(p); # is actually automatic
	class(p) = c("pm", class(p));
	return(p);
}

toMonom.pm = function(e, xsign = 1, env=NULL, verbose=TRUE) {
	if(is.null(env)) env = .GlobalEnv;
	m = data.frame(coeff=xsign);
	acc = list();
	while(TRUE) {
		if(length(e) == 1) {
			if(is.symbol(e)) {
				if(e == "-") {
					m$coeff = - m$coeff; # NO effect?
				} else if(e == "+") {
					# an extra "+"; # NO effect?
				} else {
					vn1 = as.character(e); # a variable name;
					m[, vn1] = 1;
				}
			} else if(is.language(e)) {
				pp = parse.epm(e, env=env); # another polynomial
				m  = mult.pm(pp, m);
			} else if(is.numeric(e) || is.complex(e)) {
				m[, "coeff"] = m[, "coeff"] * e;
			} else print(paste0("Error: ", e));
		} else {
			op = e[[1]];
			if(is.symbol(op)) {
				if(op == "*") {
					acc = c(acc, e[[2]]);
					e = e[[3]];
					if(is.call(e)) {
						pp = toMonom.pm(e, env=env, verbose=verbose);
						nLast = length(acc);
						pp2 = toMonom.pm(acc[[nLast]], env=env, verbose=verbose);
						pp = mult.pm(pp, pp2);
						# TODO: enforce multiplication;
						m  = mult.pm(pp, m);
						acc = acc[ - nLast];
						if(length(acc) == 0) break;
						e = acc[[length(acc)]];
					}
					next;
				}
				if(op == "(") {
					if(length(e) > 2) stop("Arrays NOT yet supported!");
					pp = parse.parenth.pm(e[[2]], env=env);
					m  = mult.pm(pp, m);
				} else if(op == "^") {
					pow = e[[3]];
					if( ! is.numeric(pow)) {
						pow = eval(substitute(pow, list(pow=pow)), envir=env);
						# pow = local(pow, envir=env);
						if( ! is.numeric(pow)) {
							warning(paste0("Power = ", pow, " is NOT numeric!"));
							pow = NA;
						}
					}
					if(is.numeric(e[[2]]) || is.complex(e[[2]])) {
						m[, "coeff"] = m[, "coeff"] * e[[2]]^pow;
					} else if(is.language(e[[2]]) && ! is.symbol(e[[2]])) {
						e = e[[2]];
						if(e[[1]] == "(") {
							pp = parse.parenth.pm(e[[2]], env=env);
							pp = pow.pm(pp, pow, debug=verbose);
							m  = mult.pm(pp, m);
						} else {
							print("Power of px!");
							pp = parse.epm(e, env=env);
							if(inherits(pp, "data.frame")) {
								pp = pow.pm(pp, pow, debug=verbose);
								m  = mult.pm(pp, m);
							} else {
								m[, "coeff"] = m[, "coeff"] * pp^pow;
							}
						}
					} else {
						vn1 = as.character(e[[2]]);
						m[, vn1] = pow;
					}
				} else if(op == "-") {
					m$coeff = - m$coeff;
					e = e[[2]]; next;
				} else if(op == "+") {
					e = e[[2]]; next;
				} else if(op == "/") {
					m[, "coeff"] = m[, "coeff"] / e[[3]];
					e = e[[2]]; next;
				} else if(is.call(e)) {
					pp = parse.epm(e, env=env); # another polynomial;
					if(is.character(pp)) {
						idv = match(pp, names(m));
						if(is.na(idv)) m[, pp] = 1
						else 
							m[, idv] = m[, idv] + 1;
					} else {
						m  = mult.pm(pp, m);
					}
					break;
				} else {
					vn1 = as.character(op); # a variable name;
					m[, vn1] = 1;
				}
			} else if(is.numeric(op) || is.complex(op)) {
				m[, "coeff"] = m[, "coeff"] * op;
			}
		}
		if(length(acc) == 0) break;
		e = acc[[length(acc)]];
		acc = head(acc, -1);
	}
	return(m);
}
parse.epm = function(e, env) {
	if(e[[1]] == "[") {
		p = eval(substitute(e, list(e=e)), env);
		return(p);
	}
	p = eval(substitute(e, list(e=e[[1]])), env);
	# p = local(e[[1]], list(e=e, env));
	pnames = names(e); len = length(pnames);
	if(len > 1) {
		# sequential substitution:
		for(i in seq(2, len)) {
			tmp = toPoly.pm(e[[i]], env=env);
			p = replace.pm(p, tmp, pnames[i]);
		}
	}
	return(p)
}
parse.parenth.pm = function(e, env) {
	p = toPoly.pm(e, env=env);
	return(p);
}


### Very Simple Parser for expressions
# - Splits expression into (text) Monomials:
parse.pm = function(e) {
	if( ! is.expression(e)) stop("Input must be an expression!")
	if( ! is.language(e[[1]])) return(NULL);
	e = e[[1]];
	e.txt = character(0);
	c.e = function(e, x.sign) {
		xi = if(nchar(x.sign) == 0) format(e[[3]]) else paste(x.sign, format(e[[3]]));
		c(e.txt, xi);
	}
	while(TRUE) {
		if(is.symbol(e[[1]])) {
			x.sign = paste0(e[[1]]);
			if(x.sign == "+") x.sign = ""
			else if(x.sign != "-") {
				e.txt = c(e.txt, format(e));
				break;
			}
		} else {
			print(e); break;
		}
		e.txt = c.e(e, x.sign);
		if(is.language(e[[2]])) e = e[[2]]
		else break;
		
	}
	return(e.txt);
}
