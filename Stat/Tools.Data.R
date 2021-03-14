########################
###
### Leonard Mada
### [the one and only]
###
### Data Tools
###
### draft v.0.1b


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





