########################
###
### Leonard Mada
### [the one and only]
###
### Graphics Tools
###
### draft v.0.1b


### Graphics Tools


####################
### Helper Functions

### Colours

### Find & Plot colours
find.col = function(name="red", start=1, max=30, bottom.mrg=8, ...) {
	is.col = grepl(name, colors());
	n.max = min(sum(is.col), start + max - 1);
	id = seq(start, n.max);
	name.col = colors()[is.col][id]
	x = rep(1, length(id)); names(x) = name.col;
	# set bottom margin
	old.par = par(mar=c(bottom.mrg,1,2,1) + 0.1)
		barplot(x, col=name.col, las=3, ...)
	par(old.par)
	invisible(name.col)
}
### Plot colours
plot.col = function(col, bottom.mrg=8, ...) {
	x = rep(1, length(col)); names(x) = names(col);
	# set bottom margin
	old.par = par(mar=c(bottom.mrg,1,2,1) + 0.1)
		barplot(x, col=col, las=3, ...)
	par(old.par)
	invisible()
}


######################

library(lattice) 

### Panel with Barcharts
barchart.adv = function(formula, data, main, col="steelblue", xlab="Response", ylab="Percent") {
	rhs = formula[[3]];
	if(rhs[[1]] == '|') rhs = rhs[[3]];
	lenFact = function(id) {
		if(length(rhs[[id]])> 1) {
			warning("Too many levels in formula!");
			rhs = rhs[[id]];
		}
		length(levels(eval(rhs[[id]], envir=data)));
	}
	dm = c(lenFact(2),lenFact(3));
	if(missing(main)) main = NA;
	barchart(formula, data=data, col=col,
		panel = function(x, ...) {
			panel.grid(h=-1, v=0)
			panel.barchart(x, ...)
		}, 
		par.settings = list(strip.background=list(col="lightgrey"),
						layout.heights=list(strip=1.45)),
		par.strip.text = list(col="black"),
		layout=dm, cex.axis=2, ylim=c(0,1), xlab=xlab, ylab=ylab,
		scales=list(tck=c(0.8,0.8), col="black", x=list(cex=1), y=list(cex=1)),
		main=main);
}


######################

### Examples

find.col()
find.col("green")
find.col("pale")

plot.col(heat.colors(30))


