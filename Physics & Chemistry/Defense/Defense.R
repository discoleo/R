

library(shape)


#########

# png(file="Defense.KineticMissile.png")

plot.base(xlim=c(-1,12), ylim=c(-1,12))

polygon(
	c(0,2,2.4, 5,5.5, 5,2.4, 2,0,0),
	c(6,6,7,   7,7.5, 8,8,   9,9,6), col="#808080")
text(1, 9.5, "Kinetic missile")

d = 0.75
x1 = 1
for(i in seq(5)) {
	polygon(6 + c(0,x1,x1,0,0), c(0,0,d,d,0) + 11 - 2*d*i, col=1)
}

x2 = 1
y0 = 11 - 2*5*d;
y1 = 11 - d;
polygon(6 + x1 + c(0,x2,x2,0,0), c(y0,y0,y1,y1,y0), col=1)


x3 = 2
for(i in seq(5)) {
	polygon(6 + x1 + x2 + c(0,x3,x3,0,0), c(0,0,d,d,0) + 11 - d - 2*d*i, col=1)
}

text(8, 11.8, "Armour plates")

### Arrows

ya = 10.8;
arrowTH(6 + c(0,x1), c(ya,ya), d=-0.4);
arrowTH(6 + c(x1,0), c(ya,ya), d= 0.4);
#
arrowTH(6 + x1 + c(0,x2), c(ya,ya), d=-0.4);
arrowTH(6 + x1 + c(x2,0), c(ya,ya), d= 0.4);
#
arrowTH(6 + x1 + x2 + c(0,x3), c(ya,ya), d=-0.4);
arrowTH(6 + x1 + x2 + c(x3,0), c(ya,ya), d= 0.4);

yatxt = 11.1
text(6.5, yatxt, "x1", col="red")
text(7.5, yatxt, "x2", col="red")
text(9, yatxt, "x3", col="red")


# dev.off()
