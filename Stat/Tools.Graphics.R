########################
###
### Leonard Mada
### [the one and only]
###
### Graphics Tools
###
### draft v.0.1h-ref


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
findSimilar.col = function(col, tol=3, start=1, max=30, type="Saturation") {
	type = pmatch(type, c("Luminosity", "Saturation", "byLTotal", "bySTotal"));
	if(is.na(type)) stop("Type not supported!");
	### Colours
	cols = colours();
	rgbAll = col2rgb(cols);
	this = as.vector(col2rgb(col));
	### Distance
	dAll = abs(rgbAll - this);
	d1 = apply(dAll, 2, sum);
	### Type
	if(type == 1 || type == 3) {
		d2 = apply(rgbAll, 2, sum);
	} else if(type == 2 || type == 4) {
		d2 = apply(rgbAll, 2, min) - apply(rgbAll, 2, max);
	}
	if(type >= 3) {
		tmp = d1; d1 = d2; d2 = tmp;
		if(type == 3) {
			d1 = abs(d1 - sum(this));
		} else {
			d1 = abs(max(this) - min(this) + d1);
		}
	}
	if(tol[1] > 1) d1 = d1 %/% tol[1];
	if(length(tol) > 2 && tol[2] > 1) d2 = d2 %/% tol[2];
	id = order(d1, d2);
	cols = colours()[id];
	if(max > 0) cols = cols[seq(start, length.out = min(max(0, length(cols) - max), max))];
	return(cols);
}

colors.ramp = function(col1="#E02032", col2="#3220D0", middle="white", alpha=NULL) {
	if(is.null(alpha))
		return(colorRampPalette(c(col1, middle, col2), space = "rgb"));
	colorRampPalette(c(col1, middle, col2), alpha=TRUE)
}

### Plot colours
plot.col = function(col, bottom.mrg=8, ...) {
	x = rep(1, length(col));
	if(is.null(names(col))) {
		names(x) = col;
	} else {
		names(x) = names(col);
	}
	# set bottom margin
	old.par = par(mar=c(bottom.mrg,1,2,1) + 0.1)
		barplot(x, col=col, las=3, ...)
	par(old.par)
	invisible()
}

image.col = function(col, bottom.mrg=8, cex.axis=1, las=2) {
	# set bottom margin
	old.par = par(mar=c(bottom.mrg,1,bottom.mrg,1) + 0.1);
	on.exit(old.par);
	len = length(col);
	if(len > 90) {
		col = col[1:90]; len = 90;
		print("Warning: only first 90 colors shown!")
	}
	#
	nr = if(len > 60) 3 else if(len > 30) 2 else 1;
	if(nr > 1) { len = len - (len %% nr); col = col[1:len]; }
	nc = len %/% nr;
	z = matrix(seq(len), nrow=nc, ncol=nr);
	if(is.null(names(col))) {
		nms = col;
	} else {
		nms = names(col);
	}
	image(1:nc, 1:nr, z=z,
		xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "", col=col);
	id = seq(1, nc);
	axis(side = 1, at = id, labels = nms[id], 
            cex.axis = cex.axis, las = las, lwd = -1);
	if(nr > 1)
	axis(side = 3, at = id, labels = nms[id + (nr-1)*nc], 
            cex.axis = cex.axis, las = las, lwd = -1);
}


######################

library(lattice) 

### Panel with Barcharts
barchart.adv = function(formula, data, main, col=c("steelblue4", "cornflowerblue"),
		xlab="Response", ylab="Percent") {
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


### Enable/Disable Stamp
setStamp = function(set=FALSE) {
	if(set) {
		DescToolsOptions("stamp", reset=TRUE);
	} else {
		DescToolsOptions("stamp" = NULL);
		print("Stamp disabled!")
	}
	invisible();
}


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

####################

####################
### Correlations ###
####################

library(DescTools)

# based on PlotCorr in package(DescTools)
# - but done correctly;
# - many more variants in package: corrplot!
filter.corr = function(m, data=NULL, cluster=TRUE, lower=TRUE, p.max=0.10, cor.min = 0, debug=TRUE) {
	# Clustering
	if (cluster == TRUE) {
		idx = order.dendrogram(as.dendrogram(
			hclust(dist(m), method = "mcquitty") ));
		m = m[idx, idx];
		if(debug) print(idx);
    }
	if(p.max > 0) {
		if(! is.null(data)) {
			# p = PairApply(data, function(x, y) cor.test(x, y)$p.value, symmetric=TRUE);
			p = apply.pair(if(cluster) data[idx] else data,
				function(x, y) cor.test(x, y)$p.value, lower=TRUE, do.all="Symmetric");
			m[p > p.max] = NA;
		} else warning("Corr: No p-values as Data is NULL!")
	}
	if (cor.min != 0) m[abs(m) < abs(cor.min)] = NA;
	if(lower) {
		m[upper.tri(m, diag=TRUE)] = 0;
	}
	return(m);
}
expand.pair = function(n, all=FALSE) {
	if(all) return(expand.grid(n, n));
	id1 = unlist(lapply(seq(n), function(id) rep(id, n + 1 - id)));
	id2 = unlist(lapply(seq(n), function(id) seq(id, n)));
	return(cbind(id1, id2));
}
apply.pair = function(x, FUN, lower=TRUE, do.all=c("None", "Symmetric", "All"), ...) {
	if (is.matrix(x)) x = as.data.frame(x);
	#
	nc = ncol(x);
	do.all = match.arg(do.all);
	isAll  = if(do.all == "All") TRUE else FALSE;
	id = expand.pair(nc, all = isAll); 
	# enable simple lists:
	f  = function(idi) FUN(x[[id[idi, 1]]], x[[id[idi, 2]]], ...);
	m  = matrix(0, nrow=nc, ncol=nc);
	rez = sapply(seq(nrow(id)), f);
	if(do.all == "All") {
		m = rez;
	} else if(lower) {
		m[lower.tri(m, diag=TRUE)] = rez;
		if(do.all == "Symmetric") {
			m[upper.tri(m)] <- t(m)[upper.tri(m)];
		}
	} else {
		# TODO:
		m[upper.tri(m, diag=TRUE)] = rez;
	}
	dimnames(m) <- list(names(x), names(x))
	return(m)
}
# - many more variants in package: corrplot!
# TODO: diagonal labels;
plot.corr = function (m, lbl = TRUE, len = 20, cex.lbl = 0.75, cex.axis = 1,
		cols = colors.ramp()(len),
		breaks = seq(-1, 1, length = len + 1), border = "grey",
		lwd = 1, args.colorlegend = NULL, xaxt = par("xaxt"), yaxt = par("yaxt"),
		las = 2, mar = c(7, 7, 3, 7), legend.offset = 0.75,
		main = "", stamp=TRUE, ...) 
{
	par.old = par(mar = mar);
	on.exit(par(par.old));
	# Plot
	m = m[nrow(m):1, ];
	image(x = 1:nrow(m), y = 1:ncol(m), xaxt = "n", yaxt = "n", 
		z = m, frame.plot = FALSE, xlab = "", ylab = "", 
		col = cols, breaks = breaks, ...);
	# Axis
	if (xaxt != "n")
		axis(side = 1, at = 1:nrow(m), labels = rownames(m), 
			cex.axis = cex.axis, las = las, lwd = -1);
	if (yaxt != "n") 
		axis(side = 2, at = 1:ncol(m), labels = colnames(m), 
			cex.axis = cex.axis, las = las, lwd = -1);
	# Legend
	if (is.list(args.colorlegend) || is.null(args.colorlegend)) {
		len.col  = length(cols);
		x.offset = legend.offset;
		args.colorlegend1 <- list(
				labels = sprintf("%.1f", seq(-1, 1, length = len.col/2 + 1)),
				x = nrow(m) + 0.5 + x.offset,
				y = ncol(m) + 0.5, width = nrow(m)/len.col, # option: fixed width?
				height = ncol(m), cols = cols, cex = 0.8);
		if ( ! is.null(args.colorlegend)) {
			args.colorlegend1[names(args.colorlegend)] <- args.colorlegend;
		}
		do.call("ColorLegend", args.colorlegend1)
	}
	# Border
	if ( ! is.na(border)) {
        usr = par("usr");
		xleft = 0.5; xright = nrow(m) + xleft;
		ybottom = 0.5; ytop = nrow(m) + ybottom;
        rect(xleft = xleft, xright = xright, ybottom = ybottom, ytop = ytop,
			lwd = lwd, border = border);
		# usr <- par("usr") # ???
		clip(xleft, xright, ybottom, ytop);
		abline(h = seq(-2, nrow(m) + 1, 1) - 0.5,
			v = seq(1, nrow(m) + 1, 1) - 0.5, col = border, lwd = lwd)
		do.call("clip", as.list(usr))
	}
	# Labels
	if( ! is.na(lbl)) {
		m = t(m)[nrow(m):1, ];
		labels.m(m, na.form = "n.s.", cex.lbl=cex.lbl);
	}
	if ( ! is.null(DescToolsOptions("stamp"))) 
		Stamp();
	if (main != "") 
		title(main = main)
}
labels.m = function(m, digits=3, ld=0, na.form=NULL, cex.lbl=0.75) {
	xc = matrix(rep(1:nrow(m), each=ncol(m)), ncol=ncol(m));
	yc = matrix(rep(ncol(m):1, ncol(m)), ncol=ncol(m));
	txt = Format(m, d=digits, ldigits = ld, na.form = na.form);
	idx = lower.tri(xc, diag=FALSE);
	text(x=xc[idx], y=yc[idx], label=txt[idx], cex=cex.lbl, xpd=TRUE);
}

### Test
m = cor(mtcars)
m = filter.corr(m, data=mtcars)
plot.corr(m)


######################

### Examples

find.col()
find.col("green")
find.col("pale")

plot.col(heat.colors(30))

plot.col(findSimilar.col("#0032D0", type="L", tol=3, start=10))

plot.col(findSimilar.col("#0032D0", type="byLT", tol=5, start=1))

plot.col(findSimilar.col("#0032D0", type="byST", tol=5, start=1))

image.col(findSimilar.col("#B83280", max=90))

