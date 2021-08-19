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

extract.vars = function(e) {
	if(is.expression(e)) e = e[[1]];
	signs = numeric(0);
	vars  = character(0);
	isNum = logical(0); # TODO
	while(length(e) > 1) {
		if(e[[1]] == "+") {
			signs = c(signs, 1);
		} else if(e[[1]] == "-") {
			signs = c(signs, -1);
		}
		len = length(e);
		if(len == 3) {
			# TODO: if(is.numeric(e[[3]])) {...}
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
	return(list(vars=rev(vars), signs=rev(signs)));
}


### Test
e = parse(text="x+y+z+2+e")
extract.vars(e)


e = parse(text="-x+y-z+2+e")
extract.vars(e)

