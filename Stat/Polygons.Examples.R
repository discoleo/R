########################
###
### Leonard Mada
### [the one and only]
###
### Polygon Process
### Examples
###
### draft v.0.1a


### "Polygon"-Process

source("Polygons.R")

#####################
#####################

### Examples

################
### Polygons ###
################

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


####################
####################

#################
### Triangles ###
#################

### Right A:
d = c(5,4,3)
p = as.triangle.dist(d)

plot.ini(range(p)*1.25)
polygon(p)


### Right B:
d = c(4,5,3)
p = as.triangle.dist(d)

plot.ini(range(p)*1.25)
polygon(p)


### Right C:
d = c(4,3,5)
p = as.triangle.dist(d)

plot.ini(range(p)*1.25)
polygon(p)


### Obtuse C:
d = c(4,3,6)
p = as.triangle.dist(d)
plot.ini(c(0, 7), c(0, 7))
polygon(p)
### Obtuse C: larger
d = c(4,3,6.9)
p = as.triangle.dist(d)
polygon(p, border="red")


### Obtuse A:
d = c(7,4,3.1)
p = as.triangle.dist(d)

plot.ini(range(p)*1.25)
polygon(p)


### Negative x: Obtuse B
d = c(4,8,5)
p = as.triangle.dist(d)
r = range(p)*1.25;
plot.ini(r, r + 3)
polygon(p)


#####################

#####################
### Circumscribed ###

### Ex 1:
a = c(3*pi/5, pi/7); r=1;
p = as.triangle.circle(a, r=r, type="r", asX=TRUE)

plot.ini(c(-r, r), asp=1)
polygon(p)
circle(r=r, mid=c(0,0), col="red")


### Ex 2:
a = c(3*pi/5, 8*pi/7); r=1;
p = as.triangle.circle(a, r=r, type="r", asX=TRUE)

plot.ini(c(-r, r), asp=1)
polygon(p)
circle(r=r, mid=c(0,0), col="red")


### Ex 3:
# 180 * (1 - 3/10 - 1/7) = 100;
a = c(3*pi/5, 2*pi/7); r=1;
p = as.triangle.circle(a, r=r, asX=TRUE)

plot.ini(c(-r, r), asp=1)
polygon(p)
circle(r=r, mid=c(0,0), col="red")


### Ex 4:
a = c(2*pi/5, 8*pi/7); r=1;
p = as.triangle.circle(a, r=r, asX=TRUE)

plot.ini(c(-r, r), asp=1)
polygon(p)
circle(r=r, mid=c(0,0), col="red")


### Random:
n = 8; r=1;
p = rtriangle.circle(n, r=r, asX=TRUE)

plot.ini(c(-r,r), asp=1)
tmp = lapply(seq(n), function(id) polygon(p[[id]], border=id));
circle(r=r, mid=c(0,0), col="red")


### Ex 2:
n = 8; r=1;
p = rtriangle.circle(n, r=r, asX=FALSE)

plot.ini(c(-r,r), asp=1)
tmp = lapply(seq(n), function(id) polygon(p[[id]], border=id));
circle(r=r, mid=c(0,0), col="red")


################

################
### Incircle ###

###
d = 6; r = 2;
p = as.triangle.incircle(d, r=r, prop=1/2)

plot.ini(range(p)*1.25, asp=1)
polygon(p)
circle(r=r, mid=c(d/2, r), col="red")


###
d = 6; r = 2.5;
p = as.triangle.incircle(d, r=r, prop=1/2)

plot.ini(range(p)*1.25, asp=1)
polygon(p)
circle(r=r, mid=c(d/2, r), col="red")


###
d = 6; r = 2; prop = 1/3
p = as.triangle.incircle(d, r=r, prop=prop)

plot.ini(range(p)*1.25, asp=1)
polygon(p)
circle(r=r, mid=c(prop*d, r), col="red")


###
d = 6; r = 2; prop = 1/4;
p = as.triangle.incircle(d, r=r, prop=prop)

plot.ini(range(p)*1.25, asp=1)
polygon(p)
circle(r=r, mid=c(prop*d, r), col="red")


###
d = 6; r = 2; prop = 4/5;
p = as.triangle.incircle(d, r=r, prop=prop)

plot.ini(range(p)*1.25, asp=1)
polygon(p)
circle(r=r, mid=c(prop*d, r), col="red")


### Almost degenerate:
d = 7; r = 3; prop = 1/4;
p = as.triangle.incircle(d, r=r, prop=prop)

plot.ini(c(-250, 50), c(0, 300), asp=1)
polygon(p)
circle(r=r, mid=c(prop*d, r), col="red")


######################
### Polygon Transforms

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

