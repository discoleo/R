########################
###
### Leonard Mada
### [the one and only]
###
### Barycenters: MNIST
###
### draft v.0.1a


### "Imagine"
# Imagine there's no for-loops,
# I wonder if you can!


# partially based on:
# https://www.r-bloggers.com/2018/01/exploring-handwritten-digit-classification-a-tidy-analysis-of-the-mnist-dataset/
# - see the respective blog post for the full details on data reshaping;
# - this implementation uses tensorized functions;


library(ggplot2)


setwd(".../ML")


### Data

x.raw = read.csv("MNIST.Sample200.csv")

x.lbl = x.raw$label
table(x.lbl)

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
		list(val=max.q, Count=cnt)
	})
	r = as.data.frame(do.call(rbind, r))
	r$val = unlist(r$val); r$Count = unlist(r$Count);
	if(round) r$val = round(r$val, 0)
	return(r)
}
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
scale.l = function(l, q.val) {
	lapply(seq(length(l)), function(id) {
		m = l[[id]] / q.val[id]
		m[m > 1] = 1
		return(m)
	})
}
toRow.m = function(m) {
	dim = dim(m)
	if(length(dim) == 2) dim = c(dim, 1)
	gr = expand.grid(1:dim[1], 1:dim[2])
	id = rep(seq(0, dim[3]-1), each=dim[1]*dim[2])
	data.frame(id=id, x=rep(gr[,1], dim[3]), y=rep(gr[,2], dim[3]), val=as.vector(m))
}
toRow.l = function(l) {
	l.rows = lapply(l, toRow.m)
	do.call(rbind, l.rows)
}
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


### Luminosity-Scaling
# [done right]
q.l = 8 # 5.01
max.q = qcount.m(x, q = 1 - q.l/(WIDTH*WIDTH), round=TRUE)
x.sc = scale.l(x, max.q$val)

image(x.sc[[2]])
image(x.sc[[52]])

plot.mean(x.sc, x.lbl, mid=0.5, title.lbl = "Average value of each pixel in 10 MNIST digits")


# Levels of Grey
# may be slow: can pre-compute toRow.l(x.sc);
ggplot(toRow.l(x.sc), aes(val)) +
	geom_histogram()


###################


### TODO:
# - various stuff;

