########################
###
### Leonard Mada
### [the one and only]
###
### Graphics Tools
###
### draft v.0.1d


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


library(DescTools)


### Mosaic Plot based on PlotMosaic in DescTools
# - code should be ideally inside PlotMosaic;
plot.mosaic = function(tbl, main=NULL, col=NULL, ...) {
	dm  = dim(tbl);
	len = length(dm)
	if(len == 3) {
		if(is.null(main)) main = attr(tbl, "dimnames")[[len]];
		# old.par = par(mfrow = c(1, dm[len]));
		layout(matrix(c(1,1, seq(2, dm[len])), ncol = 1 + dm[len], nrow = 1));
		sapply(seq(dm[len]), function(id) {
			mar = getMar(tbl[,,id]);
			# Note: accuracy problem due to y-margin!
			if(id == 1) {
				PlotMosaic(tbl[,,id], main=main[id], col=col, ...);
			} else {
				PlotMosaic(tbl[,,id], main=main[id], cols=col, ylab="", mar=c(5.1,0,mar[1] + 4,0), ...);
			}
		})
		# par(old.par);
	} else if(len == 2) {
		if(is.null(main)) main = paste0(names(attr(tbl, "dimnames"))[1:2], collapse=" ~ ");
		PlotMosaic(tbl, main=main, cols=col, ...);
	}
	invisible()
}
# helper function
getMar = function(x) {
	inches_to_lines <- (par("mar")/par("mai"))[1]
	lab.width <- max(strwidth(colnames(x), units = "inches")) * inches_to_lines;
	xmar <- lab.width + 1
	lab.width <- max(strwidth(rownames(x), units = "inches")) * inches_to_lines
	ymar <- lab.width + 1
	return(c(xmar, ymar))
}


######################

### Examples

find.col()
find.col("green")
find.col("pale")

plot.col(heat.colors(30))


