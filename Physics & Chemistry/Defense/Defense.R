
### Defence  Technologies

# A. Technological improvements to reduce the impact/effectiveness
#    of hand-fired anti-tank missiles:
# - such missiles may motivate fighters to choose urban areas as combat zones,
#   as it is easy to hide from tanks/armored vehicles;
# - the existence of effective defense strategies would reduce the willingness
#   to choose urban locations for combat;


### Defence against High-Kinetic Projectiles

# - interlacing armored plates;
# - interlacing distance less than the diameter of the penetrating projectile;
# - the energy of the kinetic missile is split as the missile enters
#   between 2 relatively-close parallel plates;
# - the outer parts of the missile may break off, dissipating some of the energy;
# - the central part will be directed against a full-depth plate;
# - such disposition of plates may be more effective than a solid block
#   of metal of the same thickness;


##############

library(shape)

### Quasi-Package BioShapes
# path = "..."
# devtools::load_all(path)


##############

kineticMissile = function() {
	# Kinetic Missile
	polygon(
		c(0,2,2.4, 5,5.5, 5,2.4, 2,0,0),
		c(6,6,7,   7,7.5, 8,8,   9,9,6), col="#808080")
	text(1, 9.5, "Kinetic missile")

	# Armoured Plates: Left (x1)
	d = 0.75
	x1 = 1
	for(i in seq(5)) {
		polygon(6 + c(0,x1,x1,0,0), c(0,0,d,d,0) + 11 - 2*d*i, col=1)
	}

	x2 = 1
	y0 = 11 - 2*5*d;
	y1 = 11 - d;
	polygon(6 + x1 + c(0,x2,x2,0,0), c(y0,y0,y1,y1,y0), col=1)

	# Armoured Plates: Right (x3)
	x3 = 2
	for(i in seq(5)) {
		polygon(6 + x1 + x2 + c(0,x3,x3,0,0), c(0,0,d,d,0) + 11 - d - 2*d*i, col=1)
	}

	text(8, 11.8, "Armour plates")

	### Arrows
	measure();
	
	### Comparison
	textComparison();
	solidBlock(c(x1, x2, x3));
}


textComparison = function() {
	yc = 2
	text(8, yc, "better than solid bock", col="red")
	text(8, yc - 0.75, "of equivalent thickness", col="red")
}
solidBlock = function(x) {
	x0 = 7; yc = 0.5; h = 2;
	w = sum(x, x[2]) / 2;
	polygon(
		c(0, w, w, 0, 0) + x0,
		c(yc, yc, yc-h, yc-h, yc), col="black")
}
measure = function() {
	# TODO: publish the function arrowMeasure();
	arrowTH = arrowMeasure;
	dT = c(-0.45, 0.4);
	
	ya = 10.8;
	arrowTH(6 + c(0,x1), c(ya,ya), d=-0.2, dT=dT);
	arrowTH(6 + c(x1,0), c(ya,ya), d= 0.2, dT=dT);
	#
	arrowTH(6 + x1 + c(0,x2), c(ya,ya), d=-0.2, dT=dT);
	arrowTH(6 + x1 + c(x2,0), c(ya,ya), d= 0.2, dT=dT);
	#
	arrowTH(6 + x1 + x2 + c(0,x3), c(ya,ya), d=-0.2, dT=dT);
	arrowTH(6 + x1 + x2 + c(x3,0), c(ya,ya), d= 0.2, dT=dT);
	
	yatxt = 11.1
	text(6.5, yatxt, "x1", col="red")
	text(7.5, yatxt, "x2", col="red")
	text(9, yatxt, "x3", col="red")
}


# png(file="Defense.KineticMissile.png", h=600, w=600)

plot.new()
plot.window(xlim=c(-1,12), ylim=c(-1,12))
kineticMissile()

# dev.off()


####################
####################

library(rgl)

### Interlaced Armour

m = matrix(c(3,0,0,0,2,0,0,0,1/2), ncol=3)
col0 = "#326432"
dy = 2

close3d()
open3d()
# Set 1
shade3d(cube3d(m, col=col0))
shade3d(translate3d(cube3d(m, col=col0), 0, 0, 2))
shade3d(translate3d(cube3d(m, col=col0), 0, 0, 4))
shade3d(translate3d(cube3d(m, col=col0), 0, 0, 6))
# Set 2
shade3d(translate3d(cube3d(m, col=col0), 0, dy, 1))
shade3d(translate3d(cube3d(m, col=col0), 0, dy, 3))
shade3d(translate3d(cube3d(m, col=col0), 0, dy, 5))


######################

######################
### Robust Bunkers ###
######################

### Bunker Walls
# - Resistant to kinetic penetrators;

### Principle:
# - Add an inhomogeneous wall that interferes with
#   the path of the missile;
# - Inhomogeneities consist of sufficiently large blocks
#   of hard rock or of other very hard materials;
# - Hard material is glued together using weak cement / concrete;

### External Wall:
x = c(0,6); y = c(0,8);
wall = grid.squareREC(x, y, r = c(2/5, 1/3), fill = "#324464A0", type = "i")

# png("Defense.Bunker.Wall.png")
plot.base(axt = NULL)
lines(rect0(x, c(0,8), lwd=3, fill = "#A0A0A024"))
lines(wall)

# Solid Rock
arrowSimple(c(-2, 2.25), c(-1, 3.25), lwd = 4, col = "#883264B2")
text(-2.5, -1.2, "Solid rock / Hard concrete", cex = 1.5, adj = c(0,1))

# Weak Cement
arrowSimple(c(3, 1.85), c(-0.5, 0.7), lwd = 4, col = "#927248B2")
text(3, -0.4, "Weak concrete", cex = 1.5, adj = c(-0.05, 0.75))

# Elastic Wall
lines(rect0(x[2] + c(0,1), y, lwd=3, fill = "#9090A0B6"))
arrowSimple(c(-0.5, x[2] + 0.5), y[2] + c(1, -0.5), lwd = 4, col = "#927248B2")
text(-0.5, y[2] + 1, "Elastic\nreinforced concrete", cex = 1.5, adj = c(0.3, -0.1))

# Fill
lines(rect0(x[2] + c(1,2), y, lwd=3, fill = "#B2B2A0A2"))
arrowSimple(x[2] + c(1.5,1.5), c(-1.5, 0.7), lwd = 4, col = "#727272B2")
text(x[2] + 1.5, -1.5, "Filler / Sand", cex = 1.5, adj = c(0.5, 1.1))

# Proper Wall
lines(rect0(x[2] + c(2,4), y, lwd=3, fill = "#8080A0B2"))
arrowSimple(x[2] + c(2.7,2.7),y[2] + c(0.5, -1), lwd = 4, col = "#6060A0B2")
text(x[2] + 2.7, y[2] + 0.5, "Proper\nWall", cex = 1.5, adj = c(0.5, -0.25))

# Missile:
arrowSimple(c(-2, 0.5), c(4,4), d.lines = c(-1,1)/8, lwd = 2)

# dev.off()

