########################
###
### Leonard Mada
### [the one and only]
###
### Formula Tools
###
### draft v.0.2b


### Tools to Process Formulas & Expressions


###############
### History ###
###############

### draft v.0.2a - v.0.2b:
# - moved Code tools to new file: Tools.Code.R;
### draft v.0.1g - v.0.1g-improve:
# - cut.code() into code blocks;
# - small improvements; [v.0.1g-improve]
### draft v.0.1f - v.0.1f-refactor:
# - extract code tokens from R code;
# - [refactored] uniform result;
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

