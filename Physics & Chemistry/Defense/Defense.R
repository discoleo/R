
### Defense  Technologies

# Technological improvements to reduce the impact/effectiveness
# of hand-fired anti-tank missiles:
# - such missiles may motivate fighters to choose urban areas as combat zones,
#   as it is easy to hide from tanks/armored vehicles;
# - the existence of effective defense strategies would reduce the willingness
#   to choose urban locations for combat;


### Defense against High-Kinetic Projectiles

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
