

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

########################

### Simplify Expressions
simplify = function(x) {
	isFrNum = function(x) {
		if(x[[1]] != "/") return(FALSE);
		if(length(x) != 3) return(FALSE);
		if(is.numeric(x[[2]]) && is.numeric(x[[3]])) return(TRUE);
		return(FALSE);
	}
	# Simplify Sum
	simplifyDif = function(x) {
		len = length(x[[2]]);
		if(len == 1 && is.numeric(x[[2]])) {
			len2 = length(x[[3]]);
			if(len2 == 1) {
				# Type: a - b
				if(is.numeric(x[[3]])) {
					x = x[[2]] - x[[3]];
				}
			} else if(len2 == 3) {
				# Type: a - b/c
				if(isFrNum(x[[3]])) {
					fr = x[[3]];
					x[[3]][[2]] = x[[2]] * fr[[3]] - fr[[2]];
					x = x[[3]];
				}
			}
		} else if(len == 3 && x[[2]][[1]] == "/") {
			tmp = x[[2]];
			div = tmp[[3]];
			# Type: a/b - c
			if(is.numeric(div) && is.numeric(tmp[[2]])) {
				tmp[[2]] = tmp[[2]] - div * x[[3]];
				x = tmp;
			}
			# TODO: 2 fractions;
			# Type: a/b - c/d
		}
		return(x)
	}
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

### Tests:

z = expression(3 - 1)[[1]]
simplify(z)

z = expression(1/3 - 1)[[1]]
simplify(z)

z = expression(2 - 1/3)[[1]]
simplify(z)

z = expression((x^(1/3))^2)[[1]]
simplify(z)
