########################
###
### Leonard Mada
### [the one and only]
###
### Formula Tools
###
### draft v.0.1c-fix


### Tools to Process Formulas & Expressions


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
# is.function(eval(e[[1]][[2]][[2]][[1]])) # but error on "x";


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
parse.simple = function(x, eol="\n") {
	len = nchar(x);
	n.comm = list(integer(0), integer(0));
	n.str  = list(integer(0), integer(0));
	is.hex = function(ch) {
		# Note: only for 1 character!
		return((ch >= "0" && ch <= "9") ||
			(ch >= "A" && ch <= "F") ||
			(ch >= "a" && ch <= "f"));
	}
	npos = 1;
	while(npos <= len) {
		s = substr(x, npos, npos);
		if(s == "#") {
			n.comm[[1]] = c(n.comm[[1]], npos);
			while(npos < len) {
				npos = npos + 1;
				if(substr(x, npos, npos) == eol) break;
			}
			n.comm[[2]] = c(n.comm[[2]], npos);
			npos = npos + 1; next;
		}
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
		npos = npos + 1;
	}
	return(list(str = n.str, comm = n.comm));
}
extract.str = function(s, npos) {
	if(length(npos[[1]]) == 0) return(character(0));
	sapply(seq(length(npos[[1]])), function(id) substr(s, npos[[1]][[id]], npos[[2]][[id]]));
}
extract.str.fun = function(fn, pkg, type=1) {
	fn = as.symbol(fn); pkg = as.symbol(pkg);
	fn = list(substitute(pkg ::: fn));
	# deparse
	s = paste0(do.call(deparse, fn), collapse="");
	npos = parse.simple(s);
	extract.str(s, npos[[type]])
}
extract.str.pkg = function(pkg, type=1, exclude.z = TRUE) {
	nms = ls(getNamespace(pkg));
	l = lapply(nms, function(fn) extract.str.fun(fn, pkg));
	if(exclude.z) {
		hasStr = sapply(l, function(s) length(s) >= 1);
		nms = nms[hasStr];
		l = l[hasStr];
	}
	names(l) = nms;
	return(l);
}

### Example

pkg = "partitions"
ls(getNamespace(pkg))

###
extract.str.pkg(pkg)

###
fn = "restrictedparts"
extract.str.fun(fn, pkg)



###
paste0(deparse(partitions::restrictedparts), collapse=" ");
