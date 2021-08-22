########################
###
### Leonard Mada
### [the one and only]
###
### Data Tools
###
### draft v.0.1c


### Tools to Process/Transform Data

###############

### Groups / Aggregates

# TODO

### Row DF
# Matrix => Row DF
toRow.df = function(m, val=rownames(m), group=colnames(m), asNumVal=TRUE, asNumGroup=TRUE) {
	if(asNumVal && (! is.numeric(val))) val = as.numeric(val);
	if(asNumGroup && (! is.numeric(group))) group = as.numeric(group);
	count.df = data.frame(
		v = val,
		count = as.vector(m),
		group = rep(group, each=nrow(m))
	)
}

### Order
order.df = function(x, decreasing=TRUE, lim=c(1E+5, 10)) {
	if( ! is.null(lim)) {
		if(nrow(x) > lim[1] || ncol(x) > lim[2]) {
			msg = "Data frame is too big! Set lim = NULL to sort big data frames."
			stop(msg);
		}
	}
	order.s = function(...) order(..., decreasing=decreasing);
	id = do.call(order.s, x);
	return(id);
}

### Duplicates
countDuplicates = function(m, onlyDuplicates=FALSE) {
	count = rep(0, nrow(m));
	xd = aggregate(count ~ ., m, length);
	if(onlyDuplicates) {
		xd = xd[xd$count > 1, , drop=FALSE];
	}
	return(xd)
}

### Encrypt IDs
encrypt = function(x, offset=0, isRandom=TRUE, DEBUG=TRUE) {
	# TODO: multiple columns in df;
	old.id = unique(x)
	len = length(old.id)
	if(DEBUG) print(len)
	
	new.id = seq(1+offset, len+offset)
	if(isRandom) {
		new.ids = sample(new.id, len)
	} else {
		new.ids = new.id
	}
	new.id = new.ids[match(x, old.id)]
	return(new.id)
}


##########################
##########################

### Formulas / Expressions

extract.vars = function(e, unique=TRUE, debug=TRUE) {
	if(is.expression(e)) e = e[[1]];
	if(e[[1]] == '~') e = e[[2]];
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
	if(unique) {
		isDuplicated = duplicated(vars);
		if(any(isDuplicated)) {
			vars = vars[ ! isDuplicated];
			signs = signs[ ! isDuplicated];
			if(debug) print("Duplicates excluded!");
		}
	}
	return(list(vars=vars, signs=signs));
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

e = parse(text="~ -x+y-z+2|+e+x-y+z")
extract.vars(e)

