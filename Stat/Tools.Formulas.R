########################
###
### Leonard Mada
### [the one and only]
###
### Formula Tools
###
### draft v.0.1a


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


