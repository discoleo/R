#####################
###
### Team Project 2021
### West University
###
### Team Project: Modeling Spread of Infections
### Team: Calea D et al;
### Supervisor: Leonard Mada / Syonic
###
### draft v0.1a

### Basic Diffusion Functions

### Other functions:
# - raster with stratified population;
# - modify raster: generate obstacles (e.g. rivers with bridges);

# - improved functions compared to the student's project;


####################

### helper Functions

# - uses plot.rs from file:
#   Percolation.R;


### Basic raster (matrix)

newRaster = function(x=80, y=x, init=c(0,5), val=0, val.init=1) {
  # init: small y offset for better compatibility with the bridges;
  m = matrix(val, nrow=y, ncol=x)
  if( ! is.null(init)) {
    mid = round(dim(m) / 2) + init[c(2,1)]; # (row, col)
    m[mid[1], mid[2]] = val.init;
  }
  invisible(m)
}

### Populations

# initializes the matrix with 30% old people 
population = function(x=80, y=x, p=0.3) {
	len = length(p)
	if(len == 1) {
		m = matrix(rbinom(x * y, 1, p), ncol = x, nrow = y)
	} else {
		psum = sum(p);
		if(round(psum, 3) < 1) { p = c(1 - psum, p); }
		else len = len - 1;
		id = seq(0, len);
		m = matrix(sample(id, x * y, TRUE, p), ncol = x, nrow = y)
	}
  invisible(m)
}

### Constraints

### 7 Bridges of Konigsberg
# n = number of bridges;
# m = world matrix;
# w = width of bridge;
# y = position of bridge:
# - middle row if y = NULL;
# dy = offset from middle;
bridge = function(n, m, w=5, dy=0, y=NULL){
	ydim = nrow(m);
	ny = length(n);
	if(is.null(y)) {
		if(ny == 1) {
			nr = ydim %/% 2 + dy;
		} else {
			nr = (seq(1, ny) * ydim) %/% (ny+1) + dy;
		}
	} else {
		nr = y; # dy is ignored;
	}
	# x-positions
	xdim = ncol(m);
	len = round(xdim / n);
	x.mid = len %/% 2;
	if(length(w) == 1) w = rep(w, ny);
	w2 = w %/% 2; x0 = 1 + x.mid - w2 - (w %% 2);
	br.x = lapply(seq(ny),
		function(id) seq(x0[id], x.mid[id] + w2[id], by=1));
	if(ny == 1) {
		br.x = br.x[[1]];
		for(nc in seq(xdim)) {
			if((nc %% len) %in% br.x) next;
			m[nr, nc] = -1;
		}
	} else {
		for(id in seq(ny)) {
			for(nc in seq(xdim)) {
				if((nc %% len[id]) %in% br.x[[id]]) next;
				m[nr[id], nc] = -1;
			}
		}
	}
	invisible(m)
}

###################
###################

### Test

### Init Raster
xdim = 200
m0 = newRaster(x=xdim)

### Ex 1:
m = bridge(5, m0)

plot.rs(m)


### Ex 2:
m = bridge(c(5, 7), m0)

plot.rs(m)


### Ex 3:
m = bridge(c(5, 7, 5), m0)

plot.rs(m)

### Diffusion:
# - will be in a different file;

