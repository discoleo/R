

source("Chem3D.R")
source("Chem3D.Tests.R")


### Examples

### Projection on Plane:

###
xyz = cbind(x = c(0,3,2), y = c(0,0,6), z = c(0,0,0))
p   = c(2,1,3)
tmp = test.proj.plane(p, xyz)
points3d(p[1], p[2], 0, col = "green", size=7, alpha = 0.5)

### Special Case:
xyz = cbind(x = c(0,3,3), y = c(0,0,5), z = c(0,0,0))
p   = c(2,1,4)
tmp = test.proj.plane(p, xyz)

###
xyz = cbind(x = c(0,5,2), y = c(0,1,6), z = c(0,0,1))
p   = c(2,1,4)
tmp = test.proj.plane(p, xyz)

###
xyz = cbind(x = c(0,5,2), y = c(0,1,6), z = c(0,0,1))
p   = c(5,2,4)
tmp = test.proj.plane(p, xyz)

### Special Case:
# Proj(p on 1,2) = p2
xyz = cbind(x = c(0,5,2), y = c(0,1,6), z = c(0,0,1))
p   = c(5,1,4)
tmp = test.proj.plane(p, xyz)


### Distance between 2 Lines

# Black: Base-segment;
# Green: create perspective;
# Red: 1st set of line segments;
# Blue: 2nd set of segments (only virtual intersection);
#   (does NOT include the point of minimal distance)


### Special: Intersect & Orthogonal
# Note: intersection points refer to L1!
test.lines.minDist.Special1()
test.lines.minDist.Special1(z = c(1/2, 3/2), col = "#4428E0", add = TRUE)


### Special: Intersect (Non-Orthogonal)
test.lines.minDist.Special2()
test.lines.minDist.Special2(z = c(1/2, 3/2), col = "#4428E0", add = TRUE)


### Special: Orthogonal
test.lines.minDist.SpecialOrtho()
test.lines.minDist.SpecialOrtho(z = c(1/2, 3/2), col = "#4428E0", add = TRUE)


### General Case:
y = 1
test.lines.minDist.General(y=y)
test.lines.minDist.General(z = c(1, -3/2), y=y, col = "#D428B0", add = TRUE)
test.lines.minDist.General(z = c(1/2, 3/2), y=y, col = "#44A0E0", add = TRUE)


###############

### Basic Maths & Transformations

### Normal
p1 = c(1,5,7); p2 = c(3,5,7);
pL = matrix(c(-4,6,-5,8,-3,8), nrow=2)
test.rotate.ortho(p1, pL)
test.rotate.ortho(p2, pL, col.point = "#D68016")


### Shift Triangle
p3 = matrix(c(1,-4,6, 5,-5,8, 7,-3,8), nrow=3)
d = 3 * seq(4)
test.shift.poly(d, p3)

###
p4 = matrix(c(1,-4,-3, 5,-5,-6, 7,-3,1), nrow=3)
p4 = rbind(p4, cbind(0.8,-0.4,0.6) %*% p4)
d = 3 * seq(4)
test.shift.poly(d, p4)



### Rotate point
p = c(2,-3,17)
pL = matrix(c(-4,6,-5,8,-3,8), nrow=2)

test.rotate.point(p, pL)
test.rotate.point(c(1,1,-10), pL)


### Arbitrary Normal

###
p = matrix(c(1,5,1,3,2,-3), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

###
pL = matrix(c(-4,6,-5,8,-3,8), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

### Special Case:
pL = matrix(c(-4,4,-1,8,-5,4), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

### Special Case:
pL = matrix(c(-4,1,-1,5,-5,1), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

###
pL = matrix(c(-4,1,-5,3,-5,2), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)


##################
##################

### Examples: Ring

###
p = matrix(c(1,0,0, 0,1,0), nrow=2, byrow=T)
m = mesh.ring(c(2, 1.6, 1.2, 1), p)
wire3d(mesh3d(vertices = t(m$V), triangles = m$M))

#############

### Examples: Cylinder

###
p = matrix(c(1,5,1,3,2,-3), nrow=2)
m = mesh.cylinder(1, p)
points3d(m$V)
lines3d(p, col = "blue")
lines3d(rbind(p[1,], p[1,] + m$N1 * 5), col = "red")
lines3d(rbind(p[1,], p[1,] + m$N2 * 5), col = "purple")

###
m = mesh.cylinder(1, p, type = "Alternating")
wire3d(mesh3d(vertices = t(m$V), triangles = m$M))
points3d(p, size = 6, col = "red")

###
m = mesh.cylinder(1, p, type = "Full")
wire3d(mesh3d(vertices = t(m$V), triangles = m$M))
points3d(p, size = 6, col = "red")

# TODO: mesh with higher resolution near boundary;


###
p = matrix(c(1,5,1,3,2,3), nrow=2)
colnames(p) = c("x", "y", "z")
test.cylinder.line3d(1, p)


### Cylinder:
### Diagonal Section

###
p  = rbind(c(1,2,3), c(1,2,5));
px = c(1,1,-3)
points3d(cylinder.section.p2(px, p, r = 0.5, type = "Ellipse"))
rotate.xz(- pi/2); rotate.xy(-pi/2)

#
points3d(cylinder.section.p2(px, p, r = 0.5, type = "Simple", nL = 65, nC = 37))

#
points3d(cylinder.section.p2(px, p, r = 0.5, nL = 65, nC = 37))


############
### Cone ###

p = rbind(c(1,1,2), c(3,4,6))
points3d(cone.vertex.alternating(3, p), col = "red", alpha = 0.75)
points3d(cone.vertex.alternating(3 - 0.1, p, phi = 0.01), alpha = 0.6)


###
points3d(cone.vertex.adaptive(3, p, nC = 47))


#############
### Torus ###

###
p = rbind(c(1,-1,3), c(4,1,1))
v = mesh.vertex.torus(p, r = 1/3, nL = 40)
points3d(v$V)
lines3d(p, col = "red")
points3d(p, col = "orange", size = 6)

###
p = p + 2 * rep(v$N2, each = 2)
v = mesh.vertex.torusAxis(p, R = -2, r = 1/2, t = 1/2, phi = c(-2*pi/3, 2*pi/3))
points3d(v$V)
lines3d(p, col = "#F08032")

#
rotate.xz(0.7)
rotate.xy(0.25)


### Common Axis
p = rbind(c(1,3,3), c(14,15,1))
v = mesh.vertex.torusAxis(p, R = 4, t = 1/4, phi = c(0, 2*pi/3))
points3d(v$V)
v = mesh.vertex.torusAxis(p, R = 4, t = 1/2, phi = c(pi/3, pi))
points3d(v$V)
v = mesh.vertex.torusAxis(p, R = 4, t = 3/4, phi = c(2*pi/3, 4*pi/3))
points3d(v$V)
lines3d(p, col = "red")


#############
### Helix ###

###
p = matrix(c(1,5,1,3,2,-3), nrow=2)
N = helix(r = 1, p) # TODO

# Simple Helix Ribbon/Strands
v = helix.ribbon(p, N, type = "Screw", dash=F)
plot.helix.ribbon(v, col = "blue")


#############

### Examples: Tetrahedron

close3d()
open3d()
Th4.base(N = 32)


###################

### Icosahedron ###

p = rbind(c(1,1,2), c(3,4,6))
v = test.icosa(p)

###
r = 1.5
v = vertex.icosa.r(r = r, N = c(1, 1, sqrt(2))/2, center = c(1,1,2))
plot.bb.icosa(v)
test.math.icosa(v)


###
close3d()
p = rbind(c(1,1,2), c(3,4,6))
v = icosa(p)
test.icosa.vset(v)

#
close3d()
v = icosa(p)
m = mesh.icosa.vset(v)
test.icosa.mesh(rm.vertex.mesh(c(1,2), m))


##################

### Torus: Cyclo-6

# Top NMR specialists have gathered in an emergency meeting
# to clarify which conformation is adopted: chair or boat?
close3d()
open3d()
test.cyclo.cylinder(19)


#####################
#####################

### Ball Bearings

# the Balls:
N = matrix(c(1,0,0, 0,1,0), ncol=3, byrow=T)
s = circles.InFixedCircle(13, 5, N)
plot.spheres.chain(s)


# Bearing:
# - Intended to run inside a tube;
r = 5
N = matrix(c(1,0,0, 0,1,0, 0,0,1), ncol=3, byrow=T)
b = bearing(r, N)
plot.bearing(b)

