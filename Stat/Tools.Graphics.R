########################
###
### Leonard Mada
### [the one and only]
###
### Data Tools
###
### draft v.0.1a


### Graphics Tools


# helper function
find.col = function(name="red", start=1, max=30, ...) {
	is.col = grepl(name, colors());
	n.max = min(sum(is.col), start + max - 1);
	id = seq(start, n.max);
	name.col = colors()[is.col][id]
	x = rep(1, length(id)); names(x) = name.col;
	barplot(x, col=name.col, las=3, ...)
}



