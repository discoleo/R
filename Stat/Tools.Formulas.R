########################
###
### Leonard Mada
### [the one and only]
###
### Formula Tools
###
### draft v.0.1f


### Tools to Process Formulas & Expressions


###############
### History ###
###############


### draft v.0.1f:
# - extract code tokens from R code;
### draft v.0.1e - v.0.1e-fix2:
# - improved version of ifelse();
# - [fixed] constant value for 1st FUN; [v.0.1e-fix2]


##########################
##########################

### Formulas / Expressions


### Extract variables & signs
extract.vars = function(e, unique=TRUE, simplify="Skip", debug=TRUE) {
	if(is.character(unique)) simplify = unique;
	simplify = pmatch(simplify, c("Skip", "First", "Last", "Aggregate"));
	if(is.na(simplify)) stop("Invalid argument simplify!")
	# "Skip" = as.unique;
	if(is.expression(e)) e = e[[1]];
	if(e[[1]] == '~') {
		e = if(length(e) == 2) e[[2]] else e[[3]]; # Note: discards e[[2]]
	}
	signs = numeric(0);
	vars  = character(0);
	isNum = logical(0); # TODO
	while(length(e) > 1) {
		if(e[[1]] == "+") {
			signs = c(signs, 1);
		} else if(e[[1]] == "-") {
			signs = c(signs, -1);
		} else if(e[[1]] == "|") {
			warning("Operator | not fully supported!");
			tmp1 = extract.vars(e[[3]], unique=unique, debug=debug);
			signs = c(signs, tmp1$signs);
			vars = c(vars, tmp1$vars);
			e = e[[2]]; next;
		}
		len = length(e);
		if(len == 3) {
			# TODO: if(is.numeric(e[[3]])) {...}
			if(is.numeric(e[[3]])) warning("Numeric values not yet implemented!")
			# as.character needed for vector vs list
			vars = c(vars, as.character(e[[3]]));
			e = e[[2]];
		} else if(len == 2) {
			vars = c(vars, as.character(e[[2]])); e = NULL; break;
		} else {
			warning("Not yet supported!"); break;
		}
	}
	if( ! is.null(e)) {
		signs = c(signs, 1);
		vars = c(vars, as.character(e));
	}
	vars = rev(vars); signs = rev(signs);
	if(simplify == 1 && unique == TRUE) simplify = 2;
	return(filter.vars(vars=vars, signs=signs, simplify=simplify, debug=debug));
}
filter.vars = function(vars, signs, simplify="First", FUN=sum, debug=TRUE) {
	if(is.character(simplify)) {
		simplify = pmatch(simplify, c("Skip", "First", "Last", "Aggregate"));
		if(is.na(simplify)) stop("Invalid argument simplify!")
	}
	filterDuplicates = function(reverse=FALSE) {
		isDuplicated = if(reverse) rev(duplicated(rev(vars))) else duplicated(vars);
		if(any(isDuplicated)) {
			vars = vars[ ! isDuplicated];
			signs = signs[ ! isDuplicated];
			if(debug) print("Duplicates excluded!");
		}
		return(list(vars=vars, signs=signs));
	}
	#
	if(simplify == 1) {
		return(list(vars=vars, signs=signs));
	}
	if(simplify == 2) {
		return(filterDuplicates());
	} else if(simplify == 3) {
		return(filterDuplicates(reverse=TRUE));
	} else {
		tmp = data.frame(vars=vars, signs=signs);
		tmp = aggregate(signs ~ vars, tmp, FUN=FUN);
		tmp = tmp[tmp$signs != 0,]
		# initial order:
		id = match(tmp$vars, vars);
		id = sort(id); id = match(vars[id], tmp$vars); tmp = tmp[id, ];
		return(tmp);
	}
}

# tokenizes a formula in its parts delimited by "~"
# Note:
# - tokenization is automatic for ",";
# - but call MUST then use FUN(expression(_conditions_), other_args, ...);
split.formula = function(e) {
	tok = list();
	while(length(e) > 0) {
		if(e[[1]] == "~") {
			if(length(e) == 2) { tok = c(NA, e[[2]], tok); break; }
			tok = c(e[[3]], tok);
			e = e[[2]];
		} else {
			tok = c(e, tok); break;
		}
	}
	return(tok);
}
### Improved: ifelse()
# - multiple conditions;
# - evaluates strictly only the values for each condition;
eval.by.formula = function(e, FUN.list, ..., default=NA) {
	tok = split.formula(e);
	if(length(tok) == 0) return();
	FUN = FUN.list;
	# Argument List
	clst = substitute(as.list(...))[-1];
	len  = length(clst);
	clst.all = lapply(clst, eval);
	eval.f = function(idCond) {
		sapply(seq(length(isEval)), function(id) {
			if(isEval[[id]] == FALSE) return(default);
			# apply FUN
			args.l = lapply(clst.all, function(a) if(length(a) == 1) a else a[[id]]);
			do.call(FUN[[idCond]], args.l);
		});
	}
	# eval 1st condition:
	isEval = eval(tok[[1]]);
	if( ! is.function(FUN[[1]])) {
		rez = rep(default, length(isEval));
		rez[isEval] = FUN[[1]];
	} else rez = eval.f(1);
	if(length(tok) == 1) return(rez);
	# eval remaining conditions
	isEvalAll = isEval;
	for(id in seq(2, length(tok))) {
		if(tok[[id]] == ".") {
			# Remaining conditions: tok == ".";
			# makes sens only on the last position
			if(id < length(tok)) warning("\".\" is not last!");
			isEval = ! isEvalAll;
			if( ! is.function(FUN[[id]])) rez[isEval] = FUN[[id]]
			else rez[isEval] = eval.f(id)[isEval];
			next;
		}
		isEval = rep(FALSE, length(isEval));
		isEval[ ! isEvalAll] = eval(tok[[id]])[ ! isEvalAll];
		isEvalAll[isEval] = isEval[isEval];
		if( ! is.function(FUN[[id]])) rez[isEval] = FUN[[id]]
		else rez[isEval] = eval.f(id)[isEval];
	}
	return(rez);
}

##################

#############
### Tests ###
#############

### Ifelse variant

x = 1:10
FUN = list(function(x, y) { x*y; }, function(x, y) { x^2; }, 0);
eval.by.formula((x > 5 & x %% 2) ~ (x <= 5) ~ ., FUN, y=2, x)
eval.by.formula((x > 5 & x %% 2) ~ (x <= 5) ~ ., FUN, x=2, x)
eval.by.formula((x > 5 & x %% 2) ~ (x <= 5) ~ ., FUN, x, y=x-1)
FUN = list(FUN[[1]], 1, 0);
eval.by.formula((x > 5 & x %% 2) ~ (x <= 5) ~ ., FUN, y=2, x)
FUN = list(0, FUN[[1]], -1);
eval.by.formula((x > 5 & x %% 2) ~ (x <= 5) ~ ., FUN, x, y=x+1)


################
### Extract Sign

### Test
e = parse(text="x+y+z+2+e")
extract.vars(e)

e = parse(text="-x+y-z+2+e")
extract.vars(e)

e = parse(text="-x+y-z+2+e+x-y")
extract.vars(e)

e = parse(text="~ -x+y-z+2+e+x-y")
extract.vars(e)

e = parse(text="~ -x+y-z+2+e+x-y")
extract.vars(e, simplify="Last")

e = parse(text="~ -x+y-z+2+e+x-y")
extract.vars(e, simplify="Aggregate")

e = parse(text="~ -x+y-z+2|+e+x-y+z")
extract.vars(e)


### Experimental: Functions
e = parse(text = "~ median(x) + (min(x) - max(x))")
is.name(e[[1]][[2]][[2]][[1]])
ef = e[[1]][[2]][[2]];
### Function
ef[[1]]
if(is.name(ef[[1]]) && length(ef) > 1) print("Function!");
### Argument:
ef[[2]]
# is.function(eval(ef[[1]]));


####################
####################

### TODO
# - Explore: remindR;

### Extract/Summary Args

summary.args = function(e) {
	rs = lapply(e, function(e) {
		if(length(e) == 0) return(data.frame(type="NULL"));
		if(length(e) == 1) {
			if(is.symbol(e)) {
				return(data.frame(type="No Default"));
			}
		}
		if(is.call(e)) {
			if(e[[1]] == "c") return(data.frame(type="Default val"));
			if(e[[1]] == "list") return(data.frame(type="Default val"));
			if(class(e) == "if") return(data.frame(type="Code"));
			return(data.frame(type="Call"));
		}
		if(is.character(e) && e == "") return(data.frame(type="Default val: Empty"));
		return(data.frame(type="Default val"));
	});
	nm = names(e);
	rs = do.call(rbind, rs);
	rs$Name = nm; rs = rs[, c(2,1)];
	rs$type[rs$Name == "..."] = "hasDot";
	rownames(rs) = NULL;
	return(rs);
}

# Arguments for all functions in a package:
summary.all.args = function(nm) {
	f = ls(getNamespace(nm))
	r = lapply(seq_along(f), function(id) {
		fn = f[id];
		# if(substr(fn,1,1) %in% c("[", "_"))
		# DONE also: "<-"
		fn = paste0("\"", fn, "\"");
		e = parse(text=paste0("formals(", nm, ":::", fn, ")"));
		e = eval(e);
		if(is.null(e)) return(data.frame(Name=NA, type=NA, FUN=f[id]));
		a = summary.args(e);
		a$FUN = f[id];
		return(a);
	})
	
	do.call(rbind, r);
}


e = formals(stats:::plot.lm)
a = summary.args(e)
aggregate(rep(1, nrow(a)) ~ type, a, FUN=length)


f = ls(getNamespace("partitions"))
e = parse(text=paste0("formals(partitions:::", f[39], ")"))
e = eval(e)
a = summary.args(e)
aggregate(rep(1, nrow(a)) ~ type, a, FUN=length)


summary.all.args("partitions")


#####################
#####################


# minimalistic parser:
parse.simple = function(x, eol="\n", all.tokens=FALSE) {
	len = nchar(x);
	n.comm = list(integer(0), integer(0));
	n.str  = list(integer(0), integer(0));
	if(all.tokens) {
		# Type: 1 = "{", 2 = "(", 3 = "[";
		tk.df = data.frame(nS=integer(0), nE=integer(0), id=integer(0),
			Type=integer(0), Nested=logical(0));
	} else tk.df = NULL;
	#
	is.hex = function(ch) {
		# Note: only for 1 character!
		return((ch >= "0" && ch <= "9") ||
			(ch >= "A" && ch <= "F") ||
			(ch >= "a" && ch <= "f"));
	}
	npos = 1;
	while(npos <= len) {
		s = substr(x, npos, npos);
		# State: COMMENT
		if(s == "#") {
			n.comm[[1]] = c(n.comm[[1]], npos);
			while(npos < len) {
				npos = npos + 1;
				if(substr(x, npos, npos) == eol) break;
			}
			n.comm[[2]] = c(n.comm[[2]], npos);
			npos = npos + 1; next;
		}
		# State: STRING
		if(s == "\"" || s == "'") {
			n.str[[1]] = c(n.str[[1]], npos);
			while(npos < len) {
				npos = npos + 1;
				se = substr(x, npos, npos);
				if(se == "\\") {
					npos = npos + 1;
					# simple escape vs Unicode:
					if(substr(x, npos, npos) != "u") next;
					len.end = min(len, npos + 4);
					npos = npos + 1;
					isAllHex = TRUE;
					while(npos <= len.end) {
						se = substr(x, npos, npos);
						if( ! is.hex(se)) { isAllHex = FALSE; break; }
						npos = npos + 1;
					}
					if(isAllHex) next;
				}
				if(se == s) break;
			}
			n.str[[2]] = c(n.str[[2]], npos);
			npos = npos + 1; next;
		}
		# brackets
		if(all.tokens) {
			type = 0;
			if(s == "{") {
				type = 1;
			} else if(s == "(") {
				type = 2;
			} else if(s == "[") {
				type = 3;
				# Note: "[[";
			}
			if(type > 0) {
				nr = nrow(tk.df);
				tk.df[nr + 1, ] = data.frame(npos, NA, nr + 1, type, FALSE);
				if(nr > 1 && is.na(tk.df$nE[nr])) tk.df$Nested[nr] = TRUE;
				npos = npos + 1; next;
			}
			#
			type = 0;
			if(s == "}") {
				type = 1;
			} else if(s == ")") {
				type = 2;
			} else if(s == "]") {
				type = 3;
				# Note: "]]";
			}
			if(type > 0) {
				id = tail(tk.df$id[is.na(tk.df$nE) & tk.df$Type == type], 1);
				tk.df$nE[id] = npos;
				npos = npos + 1; next;
			}
		}
		npos = npos + 1;
	}
	return(list(str = n.str, comm = n.comm, tokens=tk.df));
}
extract.str = function(s, npos, strip=FALSE) {
	if(length(npos[[1]]) == 0) return(character(0));
	strip.FUN = if(strip) {
			function(id) {
				if(npos[[1]][[id]] + 1 < npos[[2]][[id]]) {
					nStart = npos[[1]][[id]] + 1;
					nEnd = npos[[2]][[id]] - 1; # TODO: Error with malformed string
					return(substr(s, nStart, nEnd));
				} else {
					return("");
				}
			}
		} else function(id) substr(s, npos[[1]][[id]], npos[[2]][[id]]);
	sapply(seq(length(npos[[1]])), strip.FUN);
}
extract.str.fun = function(fn, pkg, type=1, strip=TRUE) {
	fn = as.symbol(fn); pkg = as.symbol(pkg);
	fn = list(substitute(pkg ::: fn));
	# deparse
	s = paste0(do.call(deparse, fn), collapse="");
	npos = parse.simple(s, all.tokens = (type > 2));
	extract.str(s, npos[[type]], strip=strip)
}
extract.str.pkg = function(pkg, type=1, exclude.z = TRUE, strip=TRUE) {
	nms = ls(getNamespace(pkg));
	l = lapply(nms, function(fn) extract.str.fun(fn, pkg, type=type, strip=strip));
	if(exclude.z) {
		hasStr = sapply(l, function(s) length(s) >= 1);
		nms = nms[hasStr];
		l = l[hasStr];
	}
	names(l) = nms;
	return(l);
}

###########

### Example

pkg = "partitions"
ls(getNamespace(pkg))

###
extract.str.pkg(pkg)

###
extract.str.pkg("sp") # relatively fast
# 1st run: ages! afterwards: 0.2 s;
extract.str.pkg("NetLogoR")

### All strings in "base"
extract.str.pkg("base")

###
fn = "restrictedparts"
extract.str.fun(fn, pkg)


###################

### All code tokens

### Ex 1:
extract.str.fun("compositions", "partitions", type=3)

### Ex 2:
extract.str.fun("summary.Spatial", "sp", type=3)

### Ex 3:
extract.str.fun("StrSpell", "DescTools", type=3)
extract.str.fun("StrTrim", "DescTools", type=3)

### TODO:
# - functions to filter the tokens:
#   e.g. simple calls "()", etc;
# - functions to extract the names of the called functions;


####################

### Other
paste0(deparse(partitions::restrictedparts), collapse=" ");
