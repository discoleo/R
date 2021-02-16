########################
###
### Leonard Mada
### [the one and only]
###
### Barycenters: MNIST
###
### draft v.0.2j-fix


### "Imagine"
# Imagine there's no for-loops,
# I wonder if you can!


# partially based on:
# https://www.r-bloggers.com/2018/01/exploring-handwritten-digit-classification-a-tidy-analysis-of-the-mnist-dataset/
# - see the respective blog post for the full details on data reshaping;
# - this implementation uses tensorized functions;
# - the 2nd part computes Barycenters using the Wasserstein distance;


###############

###############
### History ###

### draft v.0.2j - v.0.2j-fix:
# - more work on analysis;
# - minor fix in older code;
### draft v.0.2h - v.0.2i:
# - plot transport plan;
# - added option: mid;
### draft v.0.2g:
# - cost matrix;
# - work on / exploration of barycenter parameters;
### draft v.0.2f:
# - basic k-Means:
#   TODO: very unstable!
### draft v.0.2e:
# - plot by Group;
### draft v.0.2c - v.0.2d:
# - Normalization;
# - various subtypes of same digit;
# - more work to subclassify a digit; [v.0.2d]
### draft v.0.2a - v.0.2b-bckgr:
# - Wasserstein distance: package "Barycenter";
# - exploring package "transport";
# - exploring impact of background of the barycenter;
### draft v.0.1d - v.0.1d-tfix:
# - basic work on outliers:
#  -- extract & plot outliers;
#  -- added top() function; [v.0.1d-top]
#  -- minor fix to code; [v.0.1d-tfix]
### draft v.0.1c:
# - dist_similar.l(): dist to most similar barycenter;
### draft v.0.1b:
# - tensorized dist() and grouped tdist();
### draft v.0.1a:
# - initial functions: some tensorized base functions;


################

library(ggplot2)
library(cowplot)
library(gridGraphics)
library(patchwork)
library(magrittr)

### Wasserstein Distance
# - used in the 2nd part;
# install.packages("Barycenter")
# install.packages("transport")


setwd("C:/Users/Leo Mada/Desktop/Practica/MSc/ML")


############
### Data ###
############

x.raw = read.csv("MNIST.Sample200.csv")

x.lbl = x.raw$label
table(x.lbl)


#################
### Transform ###

### pure Data
CSV_HAS_ROWNAMES = TRUE
if(CSV_HAS_ROWNAMES) {
	x.raw = x.raw[ , -(1:2)]
} else {
	x.raw = x.raw[ , -1]
}

### Order pixels properly

order.pixels = function(x) {
	WIDTH = 28
	id = 0:(WIDTH*WIDTH - 1)
	id.new = (id %% 28) + 28*(27 - id %/% 28)
	cols = paste0("pixel", id.new)
	x.raw[ , cols]
}

x.raw = order.pixels(x.raw)

####################

### List of Matrices
WIDTH = 28
x = lapply(1:nrow(x.raw), function(id) {
		v = (unlist(x.raw[id,]))
		m = matrix(v, ncol=WIDTH, nrow=WIDTH, dimnames = NULL)
		# needed only if columns are NOT properly ordered
		# m[ , rev(seq(WIDTH))]
		m
	})

image(x[[1]])
image(x[[3]])

########################
########################

### Tensorized functions

# latest version on:
# https://github.com/discoleo/R/blob/master/Img/MNIST.Barycenters.R

### helper functions
summary.m = function(x) {
	summary(as.vector(x))
}
sum.m = function(l) {
	dim = dim(l[[1]])
	s = matrix(0, nrow=dim[1], ncol=dim[2])
	lapply(l, function(m) {
		s <<- s + m
	})
	return(s)
}
sum.grey = function(img, q=0.6) {
	f = function(img) sum(img[img < quantile(img, q)])
	if(is.list(img)) {
		sapply(img, f)
	} else {
		f(img);
	}
}
maxcount.m = function(l) {
	r = lapply(l, function(m) {
		max = max(m)
		cnt = sum(m == max)
		c(max, cnt)
	})
	return(do.call(rbind, r))
}
qcount.m = function(l, q=0.9, round=FALSE) {
	r = lapply(l, function(m) {
		max.q = quantile(m, q)
		cnt = sum(m >= max.q)
		data.frame(val=max.q, Count=cnt)
	})
	r = do.call(rbind, r)
	rownames(r) = NULL
	# r = data.table::rbindlist(r)
	if(round) r$val = round(r$val, 0)
	return(r)
}
lapply.m = function(x, FUN=median, ...) {
	# x = list of matrices;
	med = lapply(x, function(m) FUN(as.vector(m), ...))
	med = do.call(rbind, med)
	return(med)
}
match.m = function(x, m) {
	sapply(seq(nrow(m)), function(id) match(x[id], m[id,]))
}
### by group
tsum.m = function(l, group) {
	dim = dim(l[[1]])
	group = as.factor(group)
	lvl = levels(group)
	s = array(0, c(dim, length(lvl)))
	lapply(seq(length(l)), function(id) {
		s[,,group[id]] <<- s[,,group[id]] + l[[id]]
	})
	return(s)
}
dist.l = function(l, m, metric="L2") {
	dist.f = if(metric == "L2") {
		function(x.m) sqrt(sum((x.m - m)^2))
	} else if(metric == "L1") {
		function(x.m) sum(abs(x.m - m))
	} else if(metric == "Var") {
		function(x.m) sum((x.m - m)^2)
	} else {
		stop("Distance metric NOT supported!")
	}
	sapply(l, dist.f)
}
tdist.l = function(l, group, m, metric="L2", gr.offset=1) {
	# Note: group offset for digit 0;
	dist.f = if(metric == "L2") {
		function(id) sqrt(sum((l[[id]] - m[,,group[id] + gr.offset])^2))
	} else if(metric == "L1") {
		function(id) sum(abs(l[[id]] - m[,,group[id] + gr.offset]))
	} else if(metric == "Var") {
		function(id) sum((l[[id]] - m[,,group[id] + gr.offset])^2)
	} else {
		stop("Distance metric NOT supported!")
	}
	sapply(seq(length(l)), function(id) dist.f(id))
}
dist_similar.l = function(l, m, group, metric="L2", gr.offset=1, pow=1.5) {
	# Note: group offset for digit 0;
	# Note: pow used only with distance = "Lpart";
	dist.f = if(metric == "L2") {
		function(m.id, d.id) sqrt(sum((l[[d.id]] - m[,,m.id])^2))
	} else if(metric == "L1") {
		function(m.id, d.id) sum(abs(l[[d.id]] - m[,,m.id]))
	} else if(metric == "Var") {
		function(m.id, d.id) sum((l[[d.id]] - m[,,m.id])^2)
	} else if(metric == "Lpart") {
		function(m.id, d.id) sum(abs(l[[d.id]] - m[,,m.id])^pow)
	} else {
		stop("Distance metric NOT supported!")
	}
	if(missing(group)) {
		dist.all = function(id) {
			d = sapply(seq(dim(m)[3]), dist.f, id)
			d.min = min(d)
			pos = which(d == d.min)
			list(id = pos - gr.offset, d=d.min)
		}
	} else {
		dist.all = function(id) {
			d = sapply(seq(dim(m)[3]), dist.f, id)
			d.min = min(d)
			pos = which(d == d.min)
			data.frame(id = pos - gr.offset, d=d.min, dgr = d[group[id] + gr.offset])
		}
	}
	r.l = lapply(seq(length(l)), function(id) dist.all(id))
	r = do.call(rbind, r.l)
	if( ! missing(group)) {
		r$group = group
	}
	return(r)
}

### Scale
scale.l = function(l, q.val, doNormalize=FALSE) {
	l.sc = lapply(seq(length(l)), function(id) {
		m = l[[id]] / q.val[id]
		m[m > 1] = 1
		return(m)
	})
	if(doNormalize) {
		invisible(scale.norm.l(l.sc))
	} else {
		invisible(l.sc)
	}
}
scale.norm.l = function(l) {
	# TODO: "scale" also the levels of grey;
	lapply(seq(length(l)), function(id) {
		m = l[[id]] / sum(l[[id]])
		return(m)
	})
}

### Convert to Row format
toRow.m = function(m, id) {
	dim = dim(m)
	if(length(dim) == 2) dim = c(dim, 1)
	gr = expand.grid(1:dim[1], 1:dim[2])
	if(missing(id)) {
		id = seq(0, dim[3]-1);
	}
	id = rep(id, each=dim[1]*dim[2])
	data.frame(id=id, x=rep(gr[,1], dim[3]), y=rep(gr[,2], dim[3]), val=as.vector(m))
}
toRow.l = function(l, id) {
	# TODO: add automatic id?
	if(is.matrix(l)) {
		return(toRow.m(l, id))
	} else if(missing(id)) {
		l.rows = lapply(l, toRow.m)
	} else {
		l.rows = lapply(seq(length(l)), function(l.id) toRow.m(l[[l.id]], id[l.id]))
	}
	do.call(rbind, l.rows)
}
### Top samples
top = function(d, group, n) {
	top.order = order(d, decreasing=TRUE)
	id.top = tapply(top.order, group[top.order], function(id) head(id, n))
	do.call(rbind, id.top)
}
top.simple = function(d, x, n=10, decreasing=TRUE) {
	top.order = order(d, decreasing=decreasing)
	head(x[top.order], n)
}

### Plot
plot.mean = function(l, x.lbl, mid=127.5, nrow=NA, title.lbl, useTheme=TRUE) {
	### by Group
	s = tsum.m(l, x.lbl)
	### row-wise
	s.df = toRow.m(s)
	img = ggplot(data=s.df, aes(x, y, fill = val)) +
        geom_tile();
	if(is.na(nrow)) img = img + facet_wrap(~ id) else img = img + facet_wrap(~ id, nrow=nrow);
	img = img + scale_fill_gradient2(low = "white", high = "black", mid = "gray", midpoint = mid)
	if(useTheme) img = img + theme_void();
	if( ! missing(title.lbl)) img = img + labs(title = title.lbl);
	img
}
plot.mmean = function(m, m.lbl, mid=127.5, nrow=NA, title.lbl, useTheme=TRUE) {
	# m = 3D matrix with means;
	### row-wise
	s.df = toRow.m(m, id=m.lbl)
	img = ggplot(data=s.df, aes(x, y, fill = val)) +
        geom_tile();
	if(is.na(nrow)) img = img + facet_wrap(~ id) else img = img + facet_wrap(~ id, nrow=nrow);
	img = img + scale_fill_gradient2(low = "white", high = "black", mid = "gray", midpoint = mid)
	if(useTheme) img = img + theme_void();
	if( ! missing(title.lbl)) img = img + labs(title = title.lbl);
	img
}
number.groups = function(group) {
	el = sort(unique(group)); id = seq(length(el));
	tbl = table(group)
	seq.len = rep(NA, length(group))
	sapply(id, function(id) {
		isSel = group == el[id]
		seq.len[isSel] <<- seq(tbl[id])
	} )
	return(seq.len)
}
group.rows = function(x, group) {
	each.len = length(x[[1]])
	x.seq = number.groups(group)
	# inst = rep(seq(length(unique(group))), each=each.len)
	inst = as.vector(sapply(x.seq, function(id) rep(id, each.len)))
	x = toRow.l(x, id=group);
	x$inst = inst;
	invisible(x)
}
plot.group = function(x, group, mid=0.01, title="", subtitle="") {
	if(is.list(x)) {
		x = group.rows(x, group);
		group = id ~ inst;
	}
	ggplot(data=x, aes(x, y, fill = val)) +
		geom_tile(show.legend = FALSE) +
		scale_fill_gradient2(low = "white", high = "black", mid = "gray", midpoint = mid) +
		facet_grid(group) +
		labs(title = title, subtitle = subtitle) +
		theme_void() +
		theme(strip.text = element_blank())
}
plot.all = function(x, id, mid=0.01, title=NULL, subtitle=NULL) {
	# mid = 0.01 for Normalized images;
	if(missing(id)) id = if(is.list(x)) seq(length(x)) else 1;
	img = toRow.l(x, id=id)
	ggplot(data=img, aes(x, y, fill = val)) +
		geom_tile(show.legend = FALSE) +
		scale_fill_gradient2(low = "white", high = "black", mid = "gray", midpoint = mid) +
		facet_wrap(~ id) +
		labs(title = title, subtitle = subtitle) +
		theme_void() +
		theme(strip.text = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))
}

#########################
#########################

#########################
### Image Exploration ###
# [not necessary]

### Mean / Pseudo-Barycentre

### all digits
s = sum.m(x)
s / length(x)

image(s/length(x))


### by Digit
s = tsum.m(x, x.lbl)

image(s[,,1])
plot.all(lapply(seq(dim(s)[3]), function(id) s[,,id]))


### row-wise
# - plot.mean() performs automatic conversion;
s.df = toRow.m(s)
# s.df$label = rep(0:9, each=WIDTH*WIDTH)
head(s.df)

image(matrix(s.df$val[s.df$id == 1], ncol=28))

# plot directly
plot.mean(x, x.lbl, mid=127.5)


###############

### Max

max = maxcount.m(x)
head(max, 10)

# more analysis

### low max-luminosity
q.l = 6.01 # 5.01
max.q = qcount.m(x, q = 1 - q.l/(WIDTH*WIDTH), round=TRUE)
head(max.q, 20)

min(max.q$val)
max.q[max.q$Count > 100, ]
max.q[max.q$val <= 252 & max.q$Count <= 30, ]

image(x[[52]])


### BIG digits
q.l = 80 # 5.01
max.q = qcount.m(x, q = 1 - q.l/(WIDTH*WIDTH), round=TRUE)
head(max.q, 20)

x.big = top.simple(max.q[,1], x, decreasing=T)
plot.all(x.big)

x.big = top.simple(max.q[,2], x, decreasing=T)
plot.all(x.big)


### Small/Thin digits
min(max.q$val)
max.q[max.q$Count < 80, ]
max.q[max.q$val < 30, ]

image(x[[11]])

q.l = 50
min.q = qcount.m(x, q = 1 - q.l/(WIDTH*WIDTH), round=TRUE)
x.small = top.simple(min.q[,2], x, decreasing=F)
plot.all(x.small)


##################
##################

##################
### Transforms ###
##################

######################
### Luminosity-Scaling

# [done right]
q.l = 8 # 5.01
max.q = qcount.m(x, q = 1 - q.l/(WIDTH*WIDTH), round=TRUE)
# scale levels & normalize
x.sc = scale.l(x, max.q$val, doNormalize=T)

image(x.sc[[2]])
image(x.sc[[52]])

image(x[[2]] / sum(x[[2]]) - x.sc[[2]])

# png(file="Barycenters.Digits.png")
title.lbl = "Average value of each pixel in 10 MNIST digits"
plot.mean(x.sc, x.lbl, mid=0.02, title.lbl = title.lbl)

# dev.off()


### Levels of Grey
# may be slow: can pre-compute toRow.l(x.sc);
ggplot(toRow.l(x.sc), aes(val)) +
	geom_histogram()


###########################
###########################

### Minimalistic Analysis

### Distance to Barycenters
# - see major section "Barycenters" for more advanced;

### Barycenters
# - simple barycenters (average pixel value);
s.sc = tsum.m(x.sc, x.lbl)
s.sc = s.sc / length(x.sc)


### Save Barycenters
SAVE_BARYC_DATA = FALSE
if(SAVE_BARYC_DATA) {
	write.csv(s.sc, file="Barycenters.Digits.csv", row.names=FALSE)
}

### Dist
DIGIT = 1
r = dist.l(x.sc[x.lbl == DIGIT], s.sc[,,DIGIT + 1], metric="L1")
boxplot(r)


r = tdist.l(x.sc, x.lbl, s.sc, metric="L1")
boxplot(r ~ x.lbl)


### Most similar digit

### L1
r.d = dist_similar.l(x.sc, s.sc, x.lbl, metric="L1")
head(r.d, n=10)
table(r.d$id == r.d$group)

table(r.d$dgr > (r.d$d + 2))

digits.neq = r.d[r.d$dgr > (r.d$d + 2) , ]
digits.neq

# image(x.sc[[as.integer(rownames(digits.neq)[1])]])

old.par = par(mfrow=c(2, 3), mar=c(2, 2, 1, 2) + 0.1)
	sapply(1:6, function(id) image(x.sc[[as.integer(rownames(digits.neq)[id])]]) )
par(old.par)


### L...
r.d = dist_similar.l(x.sc, s.sc, x.lbl, metric="Lpart", pow=0.75)
# but is NOT better;
head(r.d, n=10)
table(r.d$id == r.d$group)

table(r.d$dgr > (r.d$d + 2))


###############

# r = tdist.l(x.sc, x.lbl, s.sc, metric="L1")
head(r)

### Outliers
isOutlier = r > 250 & x.lbl == 0;
if(any(isOutlier)) {
	image(x.sc[isOutlier][[1]])
} else {
	print("No such Outlier!")
}

### Top Outliers
top.n = 6
id.top = top(r, x.lbl, n=top.n)

head(id.top)
table(x.lbl[id.top])

# TODO: add automatic id
out.top = toRow.l(x.sc[id.top], x.lbl[id.top])
out.top$inst = rep(1:top.n, each=10*WIDTH*WIDTH)
tail(out.top)

# png(file="Barycenters.Outliers.png")
plot.group(out.top, group = inst ~ id,
	title = "Least typical digits",
	subtitle = "The 6 digits within each label that had the greatest distance to the centroid")

# dev.off()


### TODO:
# - various stuff;

######################
######################

############
### Test ###

file = "Barycenters.Digits.Big.csv"
# file = "Barycenters.Digits.csv"
x.bary = as.matrix(read.csv("Barycenters.Digits.csv"))
dim(x.bary) = c(28,28,10)

max(x.bary)

plot.mmean(x.bary, 0:9, mid=0.05, title.lbl = "Average value of each pixel in 10 MNIST digits")

### Img: Save Barycenters
SAVE_BARYC_IMG = TRUE
if(SAVE_BARYC_IMG) {
	png(file="Barycenters.Digits.Big2.png")
		plot.mmean(x.bary, 0:9, mid=0.04, title.lbl = "Average value of each pixel in 10 MNIST digits")
}
if(SAVE_BARYC_IMG) {
	dev.off()
}


#########################
#########################

###################
### Barycenters ###
###################

### Wasserstein Distance 

# Near-linear time approximation algorithms for optimal transport via Sinkhorn iteration
# https://arxiv.org/abs/1705.09634

library(Barycenter)
### using "transport"
library(transport)

cost.m = function(n, method, p=1) {
	if(length(n) > 1) n = n[1]; # TODO
	n.seq = seq(0, 1, length.out=n)
	if(missing(method)) {
		costm = as.matrix(dist(expand.grid(n.seq, rev(n.seq)), diag=TRUE, upper=TRUE))
	} else {
		costm = as.matrix(dist(expand.grid(n.seq, rev(n.seq)), method=method, diag=TRUE, upper=TRUE))
	}
	if( p != 1) costm = costm^p;
	return(costm);
}
dist.Greenkhorn = function(img1, img2, method="manhattan", p=1, ...) {
	img.dim = dim(img1);
	costm = cost.m(img.dim[2], method=method, p=p)
	len = img.dim[1] * img.dim[2];
	normalize = function(m) {s = sum(m); if(s > 1) m/s else m;}
	r1 <- matrix(img1, len, 1)
	r2 <- matrix(img2, 1, len) # assumes equal dims;
	r1 = normalize(r1); r2 = normalize(r2);
	d = Greenkhorn(r1, r2, costm=costm, ...)
	invisible(d)
}
dist.Wass = function(x, bary, FUN, ...) {
	if(missing(FUN)) {
		dist.f = function(img1, img2) {
			dist.Greenkhorn(img1, img2, ...)$Distance
		}
	} else {
		dist.f = FUN;
	}
	is.imglist = function(x) is.list(x) && is.na(match(class(bary), "pgrid"));
	if((is.imglist(bary) && length(bary) > 1) ||
		(length(dim(bary)) > 2)) {
		d = lapply(bary, function(bary) {
				cat("B ")
				dist.Wass(x, bary, FUN=dist.f, ...)
			})
		cat("\n")
		d = do.call(cbind, d)
	} else if(is.imglist(x)) {
		d = sapply(x, dist.f, img2=bary)
	} else {
		d = dist.f(x, img2=bary)
	}
	return(d)
}
background = function(x, q=0.3, type="rm") {
	# remove background
	type = match(type, c("rm", "sq", "sqq", "div"))
	if(is.na(type)) stop("Type NOT supported!")
	isBckg = x <= quantile(as.vector(x), q)
	if(type == 1) {
		x[isBckg] = 0;
	} else if(type == 2) {
		x = x^2;
	} else if(type == 3) {
		x[isBckg] = x[isBckg]^2;
	} else if(type == 4) {
		x[isBckg] = x[isBckg] / 10;
	}
	x = x / sum(x);
	return(x)
}
# multiple selection
sel = function(isDigit, type.id) {
	sel.f = function(x) {
		x[isDigit][type.id]
	}
	return(sel.f)
}
bary.f = function(x, rm.bg=FALSE, q=0.3, ...) {
	x.bary = WaBarycenter(x, ...)
	x.bary = x.bary[ , rev(seq(ncol(x.bary)))]
	# remove Background:
	# - some noise (???) in the barycenters;
	if(rm.bg) {
		x.bary = background(x.bary, q=q)
	}
	invisible(x.bary)
}
plot.tplan = function(img1, img2, tplan, plot=TRUE) {
	if(plot) {
		plot(pgrid(img1), pgrid(img2), tplan=tplan)
	} else {
		f = function() {
			old.par = par(mar=c(0,0,0,0))
			plot(pgrid(img1), pgrid(img2), tplan=tplan)
			par(old.par)
		}
		invisible(f)
	}
}
plot.alltplan = function(img1, img2, tplan, mid=0.01, useLib) {
	useLib = if(missing(useLib)) 1 else
		match(use, c("patchwork", "cowplot"));
	if(is.na(useLib)) stop("Supported Libraries: patchwork or cowplot!")
	
	if(length(mid) == 1) mid = c(mid, mid);
	img1.gg = plot.all(img1, mid=mid[1]);
	img2.gg = plot.all(img2, mid=mid[2]);
	if(missing(tplan)) tplan = transport(pgrid(img1), pgrid(img2))
	img.tpl = plot.tplan(img1, img2, tplan=tplan, plot=F)

	if(useLib == 1) {
		### using patchwork
		(img1.gg + img2.gg) / (~{img.tpl()}) + plot_layout(heights=c(1,2))
	} else {
		### using cowplot
		plot_grid(
			plot_grid(img1.gg, img2.gg, ncol=2, axis="none", labels=NULL),
			img.tpl, labels = NULL, axis="none", nrow=2, rel_heights=c(1,2))
	}
}
### IO
save.bary = function(x, digit, lambda) {
	file.name = paste0("MNIST.Barycenter.D", digit, ".L", lambda, ".csv")
	write.csv(x.bary, file=file.name, row.names=FALSE)
}

###################

# basic exploration

### Digits
DIGIT = 7;
isDigit = x.lbl == DIGIT

### Barycenters
# lambda > 50 # BUT takes long!
# lambda =  50: takes 340 s;
# lambda = 100: takes 400 s (with 24 images);
# lambda = 300: takes 400 s (but seems to be less accurate than 100);
# - Manhattan distance: helps;
lambda = 300
x.bary = bary.f(x.sc[isDigit], lambda=lambda, costm=cost.m(dim(x.sc[[1]]), method="manhattan"))

image(x.bary)

# x.bary = as.matrix(read.csv("MNIST.Barycenter.D7.L300.csv"))


# Background Noise
# - a lot of background noise in the barycenter;
summary.m(x.bary)
sum.grey(x.bary, 0.6)

### Transport Plan: using "transport"
id = 1
tr = transport(pgrid(x.sc[isDigit][[id]]), pgrid(x.bary))
plot.tplan(x.sc[isDigit][[id]], x.bary, tplan=tr)

### Noise
# - a little bit of cheating;
x2.bary = background(x.bary, type="sq")
sum.grey(x2.bary, 0.6)

id = 3
tr = transport(pgrid(x.sc[isDigit][[id]]), pgrid(x2.bary))
plot.tplan(x.sc[isDigit][[id]], x2.bary, tplan=tr)
# Q: How to "overload" local densities?
# - to move mass to locations that are already overloaded;


########
### Plot
id = 1
# png(file="MNIST.TrPlan.Ex1.png")
plot.alltplan(x.sc[isDigit][[id]], x.bary, mid=c(0.01, 0.005))

### Ex 2:
id = 3
# png(file="MNIST.TrPlan.Ex2.rmBkgr.png")
plot.alltplan(x.sc[isDigit][[id]], x2.bary)

# dev.off()


########
### Plot
# [old code]
img1 = plot.all(x.sc[isDigit][id])
img2 = plot.all(list(x2.bary))
img.tpl = plot.tplan(x.sc[isDigit][[id]], x2.bary, tplan=tr, plot=F)

### using patchwork
(img1 + img2) / (~{img.tpl()}) + plot_layout(heights=c(1,2))

### using cowplot
plot_grid(
	plot_grid(img1, img2, ncol=2, axis="none", labels=NULL),
	img.tpl, labels = NULL, axis="none", nrow=2, rel_heights=c(1,2))



SAVE_BARY = FALSE
if(SAVE_BARY) {
	save.bary(x.bary, DIGIT, lambda)
}



###############

### Wasserstein Distance
id = 1
img = x.sc[isDigit][[id]];
# img = top.img[[id]];
d = dist.Greenkhorn(img, x.bary, lambda=0.01)
d$Distance
# TODO: how to visualize ???
# image(d$Transportplan)
# - plot does NOT work with this tplan!
# plot.tplan(top.img[[id]], x.bary, tplan=d$Transportplan)


### All images of the same digit
# d = sapply(x.sc[isDigit], dist.Wass, img2=x.bary)
d = dist.Wass(x.sc[isDigit], bary=x.bary)

summary(d)

### Top outliers
top.img = top.simple(d, x.sc[isDigit])

# display top outliers
plot.all(top.img)

### Non-Digit
t1 = Sys.time();
# neg.d = sapply(x.sc[ ! isDigit], dist.Wass, img2=x.bary)
neg.d = dist.Wass(x.sc[ ! isDigit], bary=x.bary)
t2 = Sys.time();
t2 - t1;
# Time difference of 12.07101 mins;
# - for 176 images;

### Top similar
top.img = top.simple(neg.d, x.sc[ ! isDigit], n=16, decreasing=FALSE)

# display top similar images
plot.all(top.img)

### Overlap
summary(d)
summary(neg.d)

top.simple(neg.d, neg.d, n=16, decreasing=FALSE)
top.simple(neg.d, x.lbl[ ! isDigit], n=16, decreasing=FALSE)


### using "transport"
id = 1
tr = transport(pgrid(top.img[[id]]), pgrid(x.bary))
plot.tplan(top.img[[id]], x.bary, tplan=tr)


#####################
#####################

### using "transport"
library(transport)


### Top outliers
top.img = top.simple(d, x.sc[isDigit])

### Transport Outlier 1 => Barycenter
id = 1
# wasserstein(pgrid(top.img[[id]] / sum(top.img[[id]])), pgrid(x.bary))
wasserstein(pgrid(top.img[[id]]), pgrid(x.bary))
plot(pgrid(top.img[[id]]), pgrid(x.bary))

tr = transport(pgrid(top.img[[id]]), pgrid(x.bary))
plot(tr)
plot.tplan(top.img[[id]], x.bary, tplan=tr)


### Analysing the Background
# - a lot of background noise in the barycenter;
# - but NOT in the individual images;
med = lapply.m(x.sc, FUN=median)
summary(med)

med = lapply.m(x.sc, FUN=quantile, 0.8)
summary(med)

# remove background noise
x.bary = background(x.bary, 0.3)
summary(as.vector(x.bary))
image(x.bary)



##########################
##########################

### Multiple sub-classes
### of same digit

DIGIT = 7;
isDigit = x.lbl == DIGIT;

### Base-Types
### Digit 7
type.id = c(2, 17, 15)
# TODO:
# - find optimal sub-clustering "centers";

sel7 = sel(isDigit, type.id)
plot.all(sel7(x.sc), mid=0.01)

# Normalization decreases the density!
plot.all(x.sc[isDigit], mid=0.01)

### Wasserstein Distance
# - Greenkhorn: NOT accurate when using standard cost-matrix;
# d = dist.Wass(x.sc[isDigit], sel7(x.sc))
d = dist.Wass(x.sc[isDigit], sel7(x.sc),
	FUN=function(img1, img2) wasserstein(pgrid(img1), pgrid(img2)))
head(d)

### Simple
d.min = apply(d, 1, min)
head(d.min)

d.id = match.m(d.min, d)
head(d.id)
table(d.id)

### Plot by Group
plot.group(x.sc[isDigit], group=d.id)


###########
### K-Means
n.clusters = 3;
# very unstable!!!
d.cl = kmeans(d, n.clusters)
table(d.cl$cluster)

plot.group(x.sc[isDigit], group=d.cl$cluster)


###########
# TODO:

wasserstein(pgrid(sel7(x.sc)[[3]]), pgrid(sel7(x.sc)[[3]]))

x.bary = WaBarycenter(x.sc[isDigit])
x.bary = x.bary[,rev(seq(ncol(x.bary)))]


