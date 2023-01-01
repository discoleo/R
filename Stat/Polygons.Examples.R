########################
###
### Leonard Mada
### [the one and only]
###
### Polygon Process
### Examples
###
### draft v.0.1c


### "Polygon"-Process

# - various examples: Polygons;
# - Examples with Triangles have been moved
#   to file:
#   Polygons.Examples.Triangles.R;


#####################

source("Polygons.R")

#####################
#####################

################
### Polygons ###
################

### Generation:
# - from side lengths;

### Ex 1:
x = runif(10, 1, 3)
plot.dp(x, c(0.9, 0.5))

a1 = multiroot(polygonOrigin, c(0.9, 0.5), d=x)
a2 = optim(c(0.9, 0.5), polygonOptim, d=x)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a1$root)
lines.dp(x, a2$par, col="red", lty=4, lwd=2)


### Ex 2:
x = runif(20, 1, 3)
x0 = c(0.3, 0.3)
plot.dp(x, x0)

# multiroot: sensible to x0!
a1 = multiroot(polygonOrigin, x0, d=x)
a2 = optim(x0, polygonOptim, d=x)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a2$par)
lines.dp(x, a1$root, col="red", lty=4, lwd=2)


### Ex 3:
x = runif(7, 1, 4)
x0 = c(1, 0.7)
plot.dp(x, x0)

# multiroot: sensible to x0!
a1 = multiroot(polygonOrigin, x0, d=x)
a2 = optim(x0, polygonOptim, d=x)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a2$par)
lines.dp(x, a1$root, col="red", lty=4, lwd=2)


### Ex 4:
x = runif(4, 1, 5)
x0 = c(1.3, 1.3)
plot.dp(x, x0)

# multiroot: sensible to x0!
a1 = multiroot(polygonOrigin, x0, d=x)
a2 = optim(x0, polygonOptim, d=x)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a2$par)
lines.dp(x, a1$root, col="red", lty=4, lwd=2)


### Ex 5:
x = runif(4, 1, 5)
# - 3 values of x0 are necessary in some situations;
# - but only optim handles properly 3 values;
# - multiroot: very sensible to x0 with 3 values!
x0 = c(0.327, 2.78, 1.72)
plot.dp(x, x0)

a1 = multiroot(polygonOrigin, x0, d=x)
a2 = optim(x0, polygonOptim, d=x)

print(a1); cat("=========\n"); print(a2);

plot.dp(x, a2$par)
lines.dp(x, a1$root, col="red", lty=4, lwd=2)


#####################
#####################

##################
### Transforms ###
##################

### Bezier-Type

###
n = 6
t = seq(n) / (n+1)
p = transform.Bezier(t, phi=pi/5)

plot.ini(range(p)*1.25)
polygon(p)


###
n = 6
t = seq(n)*3/(4*n+1)
p = transform.Bezier(t)

plot.ini(range(p)*1.25)
polygon(p)


###
t = c(1,2,5,3,6,4)
n = length(v)
t = t / (n+1)
p = transform.Bezier(t)

plot.ini(range(p)*1.25)
polygon(p)


###
t = (1:6)^1.5
n = max(v)*2;
t = t/(n+1)
p = transform.Bezier(t)

plot.ini(range(p)*1.25)
polygon(p)
p = transform.Bezier(c(t[-1], t[1]), xy=p)
polygon(p, border="red")

