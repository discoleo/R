

### Differentiation

# Tools for Symbolic Differentiation


### Differentiate

# Inline Expression
D.expr = function(x, by) {
	y = substitute(expression(x));
	y = y[[2]];
	dx = lapply(by, function(nm) {
		D(y, nm);
	})
	return(dx)
}
# Formula
D.form = function(x, by) {
	if(length(x) == 3) {
		x = x[[3]];
	} else x = x[[2]];
	D(x, by);
}

# Partial evaluation:
eval.part = function(x, parameters) {
	do.call(substitute, list(x, parameters));
}

### Jacobian
# x = Formula;
jac.f = function(x, by) {
	ids = expand.grid(by, seq_along(x), stringsAsFactors = FALSE);
	LEN = nrow(ids);
	if(LEN == 0) return();
	lapply(seq(LEN), function(id) {
		D.form(x[[ids[id,2]]], ids[id,1]);
	})
}

########################

### Simplify Expressions
simplify = function(x) {
	# Simplify Power:
	simplifyPow = function(x) {
		if(x[[2]][[1]] == "(") {
			# TODO: LEN / ALL Tokens
			if(x[[2]][[2]][[1]] == "^") {
				tmp = x[[2]][[2]][[2]];
				pow = x[[2]][[2]][[3]];
				x[[2]] = tmp;
				if(pow[[1]] == "(") {
					if(pow[[2]][[1]] == "/") {
						num = pow[[2]][[2]];
						pow[[2]][[2]] = num * x[[3]];
						x[[3]] = pow;
					}
				} else x[[3]] = x[[3]] * pow;
			}
		}
		return(x);
	}
	n = 1;
	# TODO: traverse Tree: BFS vs DFS?
	while(TRUE) {
		if(x[[1]] == "-") {
			x = simplifyDif(x);
		} else if(x[[1]] == "^") {
			x = simplifyPow(x);
		}
		return(x)
	}
}

# Simplify Sum
simplifyDif = function(e) {
	isFrNum = function(x) {
		if(x[[1]] != "/") return(FALSE);
		if(length(x) != 3) return(FALSE);
		if(is.numeric(x[[2]]) && is.numeric(x[[3]])) return(TRUE);
		return(FALSE);
	}
	len = length(e[[2]]);
	if(len == 1) {
		if(! is.numeric(e[[2]])) return(e);
		len2 = length(e[[3]]);
		if(len2 == 1) {
			# Type: a - b
			if(is.numeric(e[[3]])) {
				e = e[[2]] - e[[3]];
			}
		} else if(len2 == 3) {
			# Type: a - b/c
			if(isFrNum(e[[3]])) {
				fr = e[[3]];
				e[[3]][[2]] = e[[2]] * fr[[3]] - fr[[2]];
				e = e[[3]];
			}
		}
	} else if(len == 3) {
		if(! isFrNum(e[[2]])) return(e);
		tmp = e[[2]];
		div = tmp[[3]];
		# Type: a/b - c
		if(is.numeric(e[[3]])) {
			tmp[[2]] = tmp[[2]] - div * e[[3]];
			e = tmp;
		} else if(isFrNum(e[[3]])) {
			# Type: a/b - c/d
			div2 = e[[3]][[3]];
			e[[1]] = as.symbol("/");
			e[[2]] = tmp[[2]]*div2 - div*e[[3]][[2]];
			e[[3]] = div * div2;
		}
	}
	return(e)
}

### Tests:

z = expression(3 - 1)[[1]]
simplify(z)

z = expression(1/3 - 1)[[1]]
simplify(z)

z = expression(2 - 1/3)[[1]]
simplify(z)

z = expression(1/2 - 1/3)[[1]]
simplify(z)
z = expression(3/4 - 2/3)[[1]]
simplify(z)

# Excluded:
z = expression(3 - x)[[1]]
simplify(z)

z = expression(3 - x/2)[[1]]
simplify(z)

z = expression(3 - 2/x)[[1]]
simplify(z)


### Pow:
z = expression((x^(1/3))^2)[[1]]
simplify(z)


###############

split.expr = function(e) {
	if(FALSE) {
	eNum = list(E = NULL, S = logical(0));
	eSym = list(E = NULL, S = logical(0));
	eFun = list(E = NULL, S = logical(0));
	} else {
		eNum = 0; eSym = 0; eFun = 0;
	}
	add.expr = function(e, add, isSum) {
		x = if(isSum) expression(1+2) else expression(1-2);
		x = x[[1]];
		x[[2]] = e;
		x[[3]] = add;
		return(x);
	}
	add.expr.new = function(e, add, isSum) {
		e$E = c(e$E, add);
		e$S = c(e$S, isSum);
		return(e);
	}
	# Init:
	e = list(e); isSumL = TRUE;
	LEN = length(e);
	while(LEN > 0) {
		tmp = e[[LEN]]; e[[LEN]] = NULL;
		isSum  = isSumL[[LEN]];
		isSumL = isSumL[-LEN]; LEN = LEN - 1;
		if(is.numeric(tmp)) {
			eNum = add.expr(eNum, tmp, isSum);
		} else if(is.symbol(tmp)) {
			eSym = add.expr(eSym, tmp, isSum);
		} else if(tmp[[1]] == "+" || tmp[[1]] == "-") {
			if(length(tmp) == 2) {
				if(tmp[[1]] == "-") isSum = ! isSum;
			}
			if(is.numeric(tmp[[2]])) {
				eNum = add.expr(eNum, tmp[[2]], isSum);
			} else if(is.symbol(tmp[[2]])) {
				eSym = add.expr(eSym, tmp[[2]], isSum);
			} else {
				tmp2 = tmp[[2]];
				if(tmp2[[1]] == "+" || tmp2[[1]] == "-") {
					e = c(e, tmp2[[2]], tmp2[[3]]);
					isSum2 = tmp2[[1]] == "+";
					isSumL = c(isSumL, isSum, isSum2);
					LEN = LEN + 2;
				} else {
					e = c(e, tmp2); isSumL = c(isSumL, isSum);
					LEN = LEN + 1;
					# print(str(tmp[[2]]));
				}
			}
			if(length(tmp) == 3) {
				e = c(e, tmp[[3]]);
				LEN = LEN + 1;
				isSum = tmp[[1]] == "+";
				isSumL = c(isSumL, isSum);
			}
		} else {
			if(is.call(tmp)) {
				eFun = add.expr(eFun, tmp, isSum);
			} else print(str(tmp));
		}
	}
	lst = list(eNum, eSym, eFun);
	return(lst);
}

### Test:
x = expression(4+2 + +3)[[1]]
split.expr(x)
x = expression(4+2 + -3)[[1]]
split.expr(x)
x = expression(4+2 - -3)[[1]]
split.expr(x)

x = expression(5+sin(x)+2+3)[[1]]
split.expr(x)

x = expression(5 + sin(x) + x + y^2 + 3*x - 5 - x/2 + pi +
	+ 2*6 + 1/5 +2+3 + sqrt(2))[[1]]
split.expr(x)

# TODO: split Calls;


