########################
###
### Leonard Mada
### [the one and only]
###
### Polymers
###
### draft v.0.1a

### Polymers

# - some experiments with random Polymers;

### Github:
# https://github.com/discoleo/R/blob/master/Stat/Polymers.R


####################
####################

### helper Functions

# - toRaster(), print.rs():
#   see file: Percolation.R;

xy0.gen = function(n, gr.dim, xy0=c(1,1)) {
	x0 = runif(n, xy0[1], gr.dim[1]);
	y0 = runif(n, xy0[2], gr.dim[2]);
	return(list(x0=x0, y0=y0, gr.dim=gr.dim));
}
rpolymer = function(n, s, alfa.range, d=5, dir0=c(0, 2*pi)) {
	n.a = s - 1;
	n.aall = n*n.a;
	a.rg.len = length(alfa.range);
	if(a.rg.len == 1) {
		alfa = sample(c(alfa.range, -alfa.range), n.aall, TRUE);
	} else {
		if(a.rg.len > 2) warning("Warning: only the first 2 values of Range will be used!")
		alfa = runif(n.aall, alfa.range[1], alfa.range[2]);
		both = rbinom(n.aall, 1, prob=c(1/2, 1/2));
		alfa[both == 1] = - alfa[both == 1];
	}
	dir0 = runif(n, dir0[1], dir0[2]);
	alfa = rbind(dir0, matrix(alfa, ncol=n));
	alfa = apply(alfa, 2, cumsum); # each col = 1 polymer;
	return(list(alfa=alfa, s=s, d=d, alfa.range=alfa.range));
}
hinges = function(pm, xy0, doConnect=TRUE) {
	s = pm$s; s1 = s + 1;
	n = ncol(pm$alfa);
	### dx
	mx = matrix(0, ncol=n, nrow=s1);
	mx[1,] = xy0$x;
	cosa = cos(pm$alfa);
	mx[-1,] = pm$d * cosa;
	mx = apply(mx, 2, cumsum);
	### dy
	my = matrix(0, ncol=n, nrow=s1);
	my[1,] = xy0$y;
	sina = sin(pm$alfa);
	my[-1,] = pm$d * sina;
	my = apply(my, 2, cumsum);
	### connect lines:
	if( ! doConnect) {
		return(invisible(list(x=mx, y=my, cosa=cosa, sina=sina)));
	}
	mxx = array(0, c(n, pm$d, s));
	for(si in seq(s)) {
		# TODO
	}
}
polymer.gen = function(pm, xy0, val=1) {
	ph = hinges(pm, xy0);
	# TODO
}
draw.polymer = function(pm, xy0, bckg=1, autoInc=TRUE) {
	ph = hinges(pm, xy0, doConnect=F);
	m = array(bckg, c(xy0$gr.dim[2], xy0$gr.dim[1], 3))
	m = as.raster(m);
	plot(m)
	len = pm$s + 1;
	n = ncol(pm$alfa);
	col = if(autoInc) {
		hex = as.hexmode(round(seq(n)*255/n));
		paste0("#00", sample(hex, n), hex);
	} else rep(0, n);
	sapply(seq(n), function(id) {
		id0 = (id-1)*len + 1; id1 = id0 + len - 1;
		lines(ph$x[id0:id1], ph$y[id0:id1], col=col[id]);
	})
	return(invisible(m));
}

#############

### Examples:

n = 40
gr.lim = c(80, 200)
xy0 = xy0.gen(n, gr.lim)

s = 10
d = 5
alfa = c(pi/3, 2*pi/3)
pm.str = rpolymer(n, s, d=d, alfa.range=alfa);
draw.polymer(pm.str, xy0);


pm = polymer.gen(pm.str, xy0);

print.rs(pm)


###################

### Example 2:
# - same x0;

n = 40
gr.lim = c(1, 200) # x0 = seq(1, 1);
xy0 = xy0.gen(n, gr.lim)

s = 10
d = 5
dir0 = c(2*pi - 2*pi/5, 2*pi + 2*pi/5)
alpha = c(pi/3, 2*pi/3)
pm.str = rpolymer(n, s, d=d, alfa.range=alfa, dir0=dir0);
#
draw.polymer(pm.str, xy0);
abline(v=1, col="red");

### Debug:
cos(pm.str$alfa[1,])
