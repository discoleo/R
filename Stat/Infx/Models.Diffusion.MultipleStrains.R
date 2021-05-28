#####################
###
### Project: Modeling Spread of Infections
###
### Team Project 2021
### West University
###
### Team: Calea D et al;
### Supervisor: Leonard Mada / Syonic
###
### L. Mada: New/Improved functionality
### draft v0.2b

### Multiple Virus Strains

####################

### helper Functions

# - uses toRaster & plot.rs from file:
#   Percolation.R;
# - specific helper functions in file:
#   Models.Diffusion.Helper.R;

sample.it = function(iter, p) {
  sample(c(iter+1, 0), 1, prob=c(p, 1-p))
}

transmit.base = function(m, iter, p, surv.time) {
	# m = matrix[nr, nc, layers]
	# layer 1: date of infection;
	# layer 2: virus type;
	m2 = m[,,1]; # copy of layer 1: infectious people;
	recov = max(1, iter - surv.time);
  
	for(nc in 1:ncol(m)) {
		for(nr in 1:nrow(m)) {
			if(m2[nr, nc] < recov) next;
			# If the subject was not previously infected with another virus
			# and a transmission event has occurred:
			# save the virus type in the 2nd layer; otherwise it stays the same.
			type = m[nr, nc, 2];
			if(nc > 1 && m[nr, nc - 1, 1] == 0) {
				isInf = sample.it(iter, p[type]);
				if(isInf > 0) {
					m[nr, nc - 1, 1] = isInf;
					m[nr, nc - 1, 2] = type;
				}
			}
			if(nc < ncol(m) && m[nr, nc + 1, 1] == 0) {
				isInf = sample.it(iter, p[type]);
				if(isInf > 0) {
					m[nr, nc + 1, 1] = isInf;
					m[nr, nc + 1, 2] = type;
				}
			}
			if(nr > 1 && m[nr - 1, nc, 1] == 0) {
				isInf = sample.it(iter, p[type]);
				if(isInf > 0) {
					m[nr - 1, nc, 1] = isInf;
					m[nr - 1, nc, 2] = type;
				}
			}
			if(nr < nrow(m) && m[nr + 1, nc, 1] == 0) {
				isInf = sample.it(iter, p[type]);
				if(isInf > 0) {
					m[nr + 1, nc, 1] = isInf;
					m[nr + 1, nc, 2] = type;
				}
			}
		}
	}
	return(m)
}


### Infection Transmission
# - mutation to V2 appears at iteration = mutate.time;
transmit = function(m, p, surv.time, mutate.time, maxIt=50, n.mutate=1) {
	for(iter in 1:maxIt) {
		if(any(mutate.time == iter)){
			# convert 1 V1 strain into V2: V1 must be still infectious;
			m = mutate.infx(m, iter, surv.time / 2, n=n.mutate);
		}
		m = transmit.base(m, iter, p, surv.time=surv.time)
	}
	invisible(m)
}
mutate.infx = function(m, iter, surv.time, n=1) {
	# convert 1 V1 strain into V2: V1 must be still infectious;
	time.inf = max(1, iter - surv.time);
	id.all = which(m[,,2] > 0 & m[,,1] >= time.inf);
	if(length(id.all) == 0) print("Sampling error!")
	else {
		id = sample(id.all, n);
		m[,,2][id] = max(m[,,2]) + 1;
	}
	invisible(m);
}


### Viruses initialization in the population matrix
newRaster.vir = function(x=80, y=x, val=0, setV=TRUE, layers=2, xy.offset=c(0,5)) {
	# Init 3d matrix with 0 (val)
	m = array(val, c(x, y, layers));
	# Init the layers with the infections starting positions
	isV = (setV == TRUE);
	if(any(isV)) {
		id = which(isV);
		len = length(id);
		print(paste0("Virus strains: ", len));
		x = (ncol(m) * seq(len)) %/% (len+1) + xy.offset[1];
		y = nrow(m) %/% 2 + xy.offset[2];
		print(x); print(y);
		m[y, x, 2] = id;
		m[y, x, 1] = 1
	}
	invisible(m)
}

###############

### Init

### Init Raster
xdim = 200;
m = newRaster.vir(x=xdim)

### Infx
surv.time = 5
# mutate.time: time when v2 starts infecting people;
p = c(0.15, 0.25, 0.25, 0.26)
mutate.time = c(100, 102, 104)

### Simulate Diffusion of Infection
system.time( {
  inf.m = transmit(m, p, surv.time=surv.time, mutate.time=mutate.time, maxIt=xdim*3) # 0.8
} )

### Plot
plot.rs(inf.m[,,1])

plot.rs(inf.m[,,1] + ifelse(inf.m[,,2] >= 2, 100, 0))
# TODO: more advanced visualisation of virus strains;

### Analysis
analyse(inf.m[,,1])
table(inf.m[,,2])


if(FALSE) {
	png(file="Epidem.Diffusion.Strains.png")
		plot.rs(inf.m[,,1] + ifelse(inf.m[,,2] == 2, 100, 0))
	dev.off()
}


#########
### Ex 2:

### Init Raster
xdim = 200;
m = newRaster.vir(x=xdim)
m[,,1] = bridge(c(5,8,5), m[,,1])
# m[,,1] = bridge(c(9,8,7,11,8,9,7), m[,,1])

### Infx
surv.time = 5
# mutate.time: time when v2 starts infecting people;
p = c(0.15, 0.25, 0.25, 0.26)
mutate.time = c(100, 102, 104)

### Simulate Diffusion of Infection
system.time( {
  inf.m = transmit(m, p, surv.time=surv.time, mutate.time=mutate.time, maxIt=xdim*3) # 0.8
} )

### Plot
plot.rs(inf.m[,,1])

plot.rs(inf.m[,,1] + ifelse(inf.m[,,2] >= 2, 100, 0))
# TODO: more advanced visualisation of virus strains;

### Analysis
analyse(inf.m[,,1])
table(inf.m[,,2])


#########
### Ex 3:

### Init Raster
xdim = 200;
m = newRaster.vir(x=xdim, setV=c(T,T,T))
# Bridges: 7 - passage for V2, but V1 partly blocked;
m[,,1] = bridge(c(5,7,5), m[,,1])
# m[,,1] = bridge(c(9,8,7,11,8,9,7), m[,,1])

### Infx
surv.time = 5
# mutate.time: time when v2 starts infecting people;
p = c(0.15, 0.11, 0.10, # initial Viruses
	0.25, 0.25, 0.26) # mutated strains
mutate.time = c(100, 102, 104)

### Simulate Diffusion of Infection
system.time( {
  inf.m = transmit(m, p, surv.time=surv.time, mutate.time=mutate.time, maxIt=xdim*3) # 0.8
} )

### Plot
plot.rs(inf.m[,,1])

### Analysis
analyse(inf.m[,,1])
table(inf.m[,,2])

