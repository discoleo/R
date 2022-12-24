


# Ideas for students


# Init
x = seq(0, 2*pi, by=0.0025);
z = cos(x) + 1i*sin(x);

######################

### Lightbulb
y = (exp(z) + exp(z^2)) / (z*exp(z^2) + exp(z))
plot(y * 1i, type="l", xlim=c(-1,1))


### Mycology
### Ascomycota & Zygomycota


### Sporangium
# - Mucor type:
#   + Columella;
# - missing Sporangiospores;
y = (exp(z) + exp(z^2)) / (z^2*exp(z^2) + exp(z))
plot(y, type="l")

# TODO:
# - the 2 "circles" actually intersect;
# - find separate equations;


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
x0 = multiroot(function(x) {
	z2 = cos(x) + 1i*sin(x);
	y = (exp(z2)*(1-z2^2) + exp(z2^2)*z2) / (exp(z2^2) - z2);
	return(Im(y));
	}, 0.87)
x0; print(x0$root, 12)


