###
### Team Project 2021
### West University
###
### Team Project: Modeling Spread of Infections
### Team: Calea D et al;
### Supervisor: Leonard Mada / Syonic
###
### draft v0.2a

### Multiple Virus Strains

####################

### helper Functions

# - uses plot.rs from file:
#   Percolation.R;

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
# - mutation to V2 appears at iteration inf.time;
transmit = function(m, p, surv.time, inf.time, maxIt=50, n.mutate=1) {
	for(iter in 1:maxIt) {
		if(iter == inf.time){
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
	id = sample(id.all, n);
	m[,,2][id] = max(m[,,2]) + 1;
	invisible(m);
}


### Viruses initialization in the population matrix
newRaster = function(x=80, y=x, val=0, setV1=TRUE, setV2=FALSE, layers=2) {
  # Init 3d matrix with 0 (val)
  m = array(val, c(x, y, layers));
  # Init the layers with the infections starting positions
  mid = round(dim(m)[1:2] / 2);
  if(setV1 && setV2) mid[1] = nrow(m) %/% 3;
  if(setV1) {
    m[mid[1], mid[2], 2] = 1
    m[mid[1], mid[2], 1] = 1
  }
  if(setV2) {
    m[nrow(m) - mid[1], mid[2], 2] = 2
    m[mid[1], mid[2], 1] = 1
  }
  invisible(m)
}

###############

### Init
surv.time = 5
# inf.time: time when v2 starts infecting people;
inf.time = 100


### Init Raster
xdim = 200;
m = newRaster(x=xdim)

### Infx
p = c(0.15, 0.25)

### Simulate Diffusion of Infection
system.time( {
  inf.m = transmit(m, p, surv.time=surv.time, inf.time=inf.time, maxIt=xdim*3) # 0.8
} )

### Plot
plot.rs(inf.m[,,1])

plot.rs(inf.m[,,1] + ifelse(inf.m[,,2] == 2, 100, 0))
# TODO: more advanced visualisation of virus strains;

if(FALSE) {
	png(file="Epidem.Diffusion.Strains.png")
		plot.rs(inf.m[,,1] + ifelse(inf.m[,,2] == 2, 100, 0))
	dev.off()
}
