########################
###
### Leonard Mada
### [the one and only]
###
### Polygon Process
### Examples: Triangles
###
### draft v.0.1a


### "Polygon"-Process

# - Examples: triangles;

#####################

source("Polygons.R")

#####################
#####################

#################
### Triangles ###
#################

### Generation:
# - from side lengths;

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

### Generation:
# - inside circumscribed circle;

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


##############
### Random ###

###
n = 8; r=1;
p = rtriangle.circle(n, r=r, asX=TRUE)

plot.ini(c(-r,r), asp=1)
circle(r=r, mid=c(0,0), col="red")
tmp = lapply(seq(n), function(id) polygon(p[[id]], border=id));


### Ex 2:
n = 8; r=1;
p = rtriangle.circle(n, r=r, asX=FALSE)

plot.ini(c(-r,r), asp=1)
circle(r=r, mid=c(0,0), col="red")
tmp = lapply(seq(n), function(id) polygon(p[[id]], border=id));


### Ex 3:
n = 8; r=1;
p = rtriangle.circle(n, r=r, type="half", asX=TRUE)

plot.ini(c(-r,r), asp=1)
circle(r=r, mid=c(0,0), col="red")
tmp = lapply(seq(n), function(id) polygon(p[[id]], border=id));


### Ex 4:
n = 8; r=1;
p = rtriangle.circle(n, r=r, type="phalf", asX=TRUE)

plot.ini(c(-r,r), asp=1)
circle(r=r, mid=c(0,0), col="red")
tmp = lapply(seq(n), function(id) polygon(p[[id]], border=id));


### Ex 5:
n = 8; r=1;
p = rtriangle.circle(n, r=r, type="random", asX=TRUE)

plot.ini(c(-r,r), asp=1)
circle(r=r, mid=c(0,0), col="red")
tmp = lapply(seq(n), function(id) polygon(p[[id]], border=id));


### Ex 6: Eq. Partitions
n = 8; r=1;
p = rtriangle.circle(n, r=r, type="eqpart", asX=TRUE)

plot.ini(c(-r,r), asp=1)
circle(r=r, mid=c(0,0), col="red")
tmp = lapply(seq(n), function(id) polygon(p[[id]], border=id));


### Ex 7: Eq. Partitions
n = 16; r=1;
p = rtriangle.circle(n, r=r, type="eqpart", asX=TRUE)

plot.ini(c(-r,r), asp=1)
circle(r=r, mid=c(0,0), col="red")
tmp = lapply(seq(n), function(id) polygon(p[[id]], border=id));


################

################
### Incircle ###

### Generation:
# - around incircle;

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

