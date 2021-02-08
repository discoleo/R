########################
###
### Leonard Mada
### [the one and only]
###
### Barycenters: MNIST
###
### draft v.0.2b-bckg


### "Imagine"
# Imagine there's no for-loops,
# I wonder if you can!


# partially based on:
# https://www.r-bloggers.com/2018/01/exploring-handwritten-digit-classification-a-tidy-analysis-of-the-mnist-dataset/
# - see the respective blog post for the full details on data reshaping;
# - this implementation uses tensorized functions;


###############

###############
### History ###

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
library(magrittr)

### Wasserstein Distance
# - used in the 2nd part;
# install.packages("Barycenter")
# install.packages("transport")


setwd(".../ML")


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
sum.m = function(l) {
	dim = dim(l[[1]])
	s = matrix(0, nrow=dim[1], ncol=dim[2])
	lapply(l, function(m) {
		s <<- s + m
	})
	return(s)
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
scale.l = function(l, q.val) {
	lapply(seq(length(l)), function(id) {
		m = l[[id]] / q.val[id]
		m[m > 1] = 1
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
	if(missing(id)) {
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
plot.group = function(x, group, mid=0.5, title="", subtitle="") {
	ggplot(data=x, aes(x, y, fill = val)) +
		geom_tile(show.legend = FALSE) +
		scale_fill_gradient2(low = "white", high = "black", mid = "gray", midpoint = mid) +
		facet_grid(group) +
		labs(title = title, subtitle = subtitle) +
		theme_void() +
		theme(strip.text = element_blank())
}
plot.all = function(x, id, mid=0.5, title="", subtitle="") {
	if(missing(id)) id=seq(length(x))
	img = toRow.l(x, id=id)
	ggplot(data=img, aes(x, y, fill = val)) +
		geom_tile(show.legend = FALSE) +
		scale_fill_gradient2(low = "white", high = "black", mid = "gray", midpoint = mid) +
		facet_wrap(~ id) +
		labs(title = title, subtitle = subtitle) +
		theme_void() +
		theme(strip.text = element_blank())
}

#####################

### Image Exploration

### Mean / Pseudo-Barycentre

### all digits
s = sum.m(x)
s / length(x)

image(s/length(x))


### by Digit
s = tsum.m(x, x.lbl)

image(s[1,,])

### row-wise
s.df = toRow.m(s)
# s.df$label = rep(0:9, each=WIDTH*WIDTH)
head(s.df)

image(matrix(s.df$val[s.df$id == 1], ncol=28))


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

min(max.q$val)
max.q[max.q$Count < 80, ]
max.q[max.q$val < 30, ]

image(x[[11]])

######################

######################
### Luminosity-Scaling
# [done right]
q.l = 8 # 5.01
max.q = qcount.m(x, q = 1 - q.l/(WIDTH*WIDTH), round=TRUE)
x.sc = scale.l(x, max.q$val)

image(x.sc[[2]])
image(x.sc[[52]])

# png(file="Barycenters.Digits.png")
plot.mean(x.sc, x.lbl, mid=0.5, title.lbl = "Average value of each pixel in 10 MNIST digits")

# dev.off()


### Levels of Grey
# may be slow: can pre-compute toRow.l(x.sc);
ggplot(toRow.l(x.sc), aes(val)) +
	geom_histogram()


###########################
###########################

### Distance to Barycenters

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

library(Barycenter)

dist.Greenkhorn = function(img1, img2) {
	img.dim = dim(img1);
	n = seq(0, 1, length.out = img.dim[2])
	costm = as.matrix(dist(expand.grid(n, rev(n)), diag=TRUE, upper=TRUE))
	len = img.dim[1] * img.dim[2];
	r1 <- matrix(img1, len, 1)
	s = sum(r1);
	if(s > 1) r1 = r1 / s;
	r2 <- matrix(img2, 1, len) # assumes equal dims;
	d = Greenkhorn(r1, r2, costm=costm)
	invisible(d)
}
dist.Wass = function(img1, img2) {
	dist.Greenkhorn(img1, img2)$Distance
}
background = function(x, q=0.3, remove=TRUE) {
	# remove background
	isBckg = x <= quantile(as.vector(x), q)
	if(remove) {
		x[isBckg] = 0
	} else {
		x[isBckg] = x[isBckg] / 10;
	}
	x = x / sum(x)
	return(x)
}

####

DIGIT = 7;
isDigit = x.lbl == DIGIT
x.bary = WaBarycenter(x.sc[isDigit])
x.bary = x.bary[,rev(seq(ncol(x.bary)))]
# remove Background
# x.bary = background(x.bary, q=0.3)

image(x.bary)

### Wasserstein Distance
d = dist.Greenkhorn(x.sc[isDigit][[1]], x.bary)
d$Distance
# TODO: how to visualize ???
image(d$Transportplan)

### All images of the same digit
d = sapply(x.sc[isDigit], dist.Wass, img2=x.bary)

summary(d)

### Top outliers
top.img = top.simple(d, x.sc[isDigit])

# display top outliers
plot.all(top.img)

### Non-Digit
t1 = Sys.time();
neg.d = sapply(x.sc[ ! isDigit], dist.Wass, img2=x.bary)
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


##############

### using "transport"
library(transport)


### Top outliers
top.img = top.simple(d, x.sc[isDigit])

### Transport Outlier 1 => Barycenter
id = 1
wasserstein(pgrid(top.img[[id]] / sum(top.img[[id]])), pgrid(x.bary))
plot(pgrid(top.img[[id]] / sum(top.img[[id]])), pgrid(x.bary))

tr = transport(pgrid(top.img[[id]] / sum(top.img[[id]])), pgrid(x.bary))
plot(tr)
plot(pgrid(top.img[[id]] / sum(top.img[[id]])), pgrid(x.bary), tplan=tr)


### Analysing the Background
# - a lot of background noise in the barycenter;
med = lapply.m(x.sc, FUN=median)
summary(med)

med = lapply.m(x.sc, FUN=quantile, 0.8)
summary(med)

summary(as.vector(x.bary))
x.bary = background(x.bary, 0.3)
image(x.bary)




