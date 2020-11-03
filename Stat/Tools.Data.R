########################
###
### Leonard Mada
### [the one and only]
###
### Data Tools
###
### draft v.0.1a


### Tools to Process/Transform Data

###############

### Encrypt IDs
encrypt = function(x, isRandom=TRUE, DEBUG=TRUE) {
	# TODO: multiple columns in df;
	old.id = unique(x)
	len = length(old.id)
	if(DEBUG) print(len)
	
	if(isRandom) {
		new.ids = sample(seq(len), len)
	} else {
		new.ids = seq(len)
	}
	new.id = new.ids[match(x, old.id)]
	return(new.id)
}





