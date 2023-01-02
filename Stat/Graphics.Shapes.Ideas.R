########################
###
### Leonard Mada
### [the one and only]
###
### Object Shapes
###
### draft v.0.1c


### Ideas for students


##################

library(rootSolve)
# used to compute some bounds;


seq.unit = function(x0=0, by=0.0025, symmetric=TRUE) {
	xEnd = if(symmetric) 2*pi - x0 else x0[2];
	x2 = seq(x0, xEnd, by=by);
	z2 = cos(x2) + 1i*sin(x2);
	return(z2);
}

# Init
x = seq(0, 2*pi, by=0.0025);
z = cos(x) + 1i*sin(x);

######################

### Light-bulb
y = (exp(z) + exp(z^2)) / (z*exp(z^2) + exp(z))
plot(y * 1i, type="l", xlim=c(-2,2), ylim=c(-2,2))
#
for(k in seq(1.125, 1.5, by=0.125)) {
	y = (exp(z) + exp(z^2)) / (z*exp(z^2) + exp(z*k));
	y = sqrt(k) * Re(y) + 1i*Im(y);
	lines(y * 1i, type="l", col="red")
}

### Chemistry: various Flasks
x0 = 1.028713;
z2 = seq.unit(x0)
flask = function(z, k=1.25, mult=1i) {
	y = (exp(z) + exp(z^2)) / (z*exp(z^2) + exp(z/k));
	y = 1/k * Re(y) + 1i*Im(y);
	return(y * mult);
}
### Flask
plot(flask(z2), type="l", xlim=c(-2,2), ylim=c(-2,2))
for(k in seq(1,1.5, by=0.125)) lines(flask(z2, k=k), col="red")
lines(flask(z2))

### Anatomy: Bladder
x0 = 1.193736;
z2 = seq.unit(x0)
plot(flask(z2, k=1.75, mult=-1i), type="l", xlim=c(-2,2), ylim=c(-2,2))


solve.x0 = function(x0, ...) {
	multiroot(function(x, k, R=1, mult=1i) {
		z = cos(x) + 1i*sin(x);
		y = (exp(z) + exp(z^2)) / (z*exp(z^2) + exp(z/k));
		y = 1/k * Re(y) + 1i*Im(y);
		y = y * mult;
		return(Im(y) - R);
	}, x0, ...)
}
#
x0 = solve.x0(1, k=1.25, R=1, mult=1i)
x0; print(x0$root, 12)
#
x0 = solve.x0(1, k=1.75, R=0, mult=1)
x0; print(x0$root, 12)


################

################
### Mycology ###
################

### Ascomycota & Zygomycota

### Sporangium
# - Mucor type:
#   + Columella;
# - missing Sporangiospores;
y = (exp(z) + exp(z^2)) / (z^2*exp(z^2) + exp(z))
plot(y * -1i, type="l", xlim=c(-2,2))

# TODO:
# - the 2 "circles" actually intersect;
# - find separate equations;

###
k = 2; # parameter
y = (exp(z) + exp(z^2/2)) / (z*exp(z) - sin(z^2/k))
plot(y * 1i, type="l", ylim=c(-4,2))
for(k in seq(2.5, 5, by=0.5)) {
	y = (exp(z) + exp(z^2/2)) / (z*exp(z) - sin(z^2/k));
	lines(y * 1i, type="l", col="red");
}


### Sporangium or Embryology
# Embryology: Transverse section;
y = (exp(z) + exp(z^2)) / (z*exp(z) - sin(1/z))
plot(y * 1i, type="l")
#
plot(y * -1i, type="l")
# Parametric
plot(y * 1i, type="l", xlim=c(-8,8), ylim=c(-4,12))
for(k in seq(1.125, 1.5, by=0.125)) {
	y2 = (exp(z) + exp(z^2)) / (z*exp(z) - sin(1/z/k));
	lines(y2 * 1i, type="l", col="red");
}


### Other
y = (exp(z) + exp(z^2)) / (sqrt(z)*exp(z^2) + exp(z))
plot(y * -1i, type="l")


### Kidney
len = length(z)
lim = round(len * 1.25/4);
z2 = z[c(seq(len - lim, len, by=1), seq(1, lim, by=1))]
y = (exp(z2) + exp(z2^2)) / (sqrt(z2)*exp(2*z2^2) + exp(z2));
y = Re(y) + 1.25i*Im(y);
plot(y, type="l", xlim=c(0.25, 1.75))


####################
####################

####################
### Parasitology ###
####################

### Parasite Egg
# - or Sporangium;
x0 = 0.87327656832
x2 = seq(x0, 2*pi - x0, by=0.00125)
z2 = cos(x2) + 1i*sin(x2)
y = (exp(z2)*(1-z2^2) + exp(z2^2)*z2) / (exp(z2^2) - z2)
plot(y, type="l", ylim=c(-8,8))

###

# Im(y) == 0;
x0 = multiroot(function(x) {
	z2 = cos(x) + 1i*sin(x);
	y = (exp(z2)*(1-z2^2) + exp(z2^2)*z2) / (exp(z2^2) - z2);
	return(Im(y));
	}, 0.87)
x0; print(x0$root, 12)




####################

### 2 Shells:
y = (exp(z) + exp(z^2/2)) / (z*exp(z) - sin(z^2/1.5))
plot(y * 1i * Im(y), type="l", ylim=c(-3,3))
#
for(k in c(1.6, 1.7, 1.8)) {
	y = (exp(z) + exp(z^2/2)) / (z*exp(z) - sin(z^2/k))
	lines(y * 1i * Im(y), type="l", col="red")
}


# Half:
k = 1.6;
x0 = pi;
x2 = seq(0, x0, by=0.0025);
z2 = cos(x2) + 1i*sin(x2);
y = (exp(z2) + exp(z2^2/2)) / (z2*exp(z2) - sin(z2^2/k))
plot(y * 1i * Im(y), type="l", ylim=c(-3,3))


# Im(y) == 0;
# - the other non-trivial "roots";
x0 = multiroot(function(x, k=1.6) {
	z = cos(x) + 1i*sin(x);
	y = (exp(z) + exp(z^2/2)) / (z*exp(z) - sin(z^2/k));
	y = y * 1i * Im(y);
	return(Im(y));
	}, 1.9) # x0 = 2*pi - 1.9;
x0; print(x0$root, 12)

