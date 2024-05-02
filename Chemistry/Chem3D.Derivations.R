

### Rotation by pi/2 around Line

# p1, p2 = points defining line;
# p0 = point which will be rotated by pi/2;
# pP = projected point on line;
# pT = rotated point;
# pR = point between (p1, p2) such that ||pR - pP|| = ||p0 - pP||;
# (may be useful in other computations)

# Note: Fraction for pR:
# pR = tR*pP + (1-tR)*p1;
tR = 1 - R / sqrt((x1 - xP)^2 + (y1 - yP)^2 + (z1 - zP)^2);


### Equations:
R^2 - (x0 - xP)^2 + (y0 - yP)^2 + (z0 - zP)^2 # = 0; # R = computable;
(x0 - xP)*(xT - xP) + (y0 - yP)*(yT - yP) + (z0 - zP)*(zT - zP) # = 0
(x1 - xP)*(xT - xP) + (y1 - yP)*(yT - yP) + (z1 - zP)*(zT - zP) # = 0
(xT - xP)^2 + (yT - yP)^2 + (zT - zP)^2 - R^2 # = 0

# Solve for (xT, yT, zT);
# dx0 = x0 - xP; dx1 = x1 - xP; dxT = xT - xP;
dx0*dxT + dy0*dyT + dz0*dzT # = 0
dx1*dxT + dy1*dyT + dz1*dzT # = 0
dxT^2 + dyT^2 + dzT^2 - R^2 # = 0
# =>
p0 = as.pm("dx0*dxT + dy0*dyT + dz0*dzT")
p1 = as.pm("dx1*dxT + dy1*dyT + dz1*dzT")
p2 = as.pm("dxT^2 + dyT^2 + dzT^2 - R^2")
pR = solve.lpm(p0,p1,p2, xn = c("dzT", "dyT"))
print.pm(pR[[2]]$Rez, leading = "dxT")

### Eq:
((dx0*dy1 - dx1*dy0)^2 + (dx0*dz1 - dx1*dz0)^2 + (dy0*dz1 - dy1*dz0)^2) * dxT^2 +
 - (dy0*dz1 - dy1*dz0)^2*R^2 # = 0

### Solutions:
dyT = - (dx0*dz1 - dx1*dz0)*dxT / (dy0*dz1 - dy1*dz0);
dzT = - (dyT*dy0 + dxT*dx0) / dz0;
dzT =   (dx0*dy1 - dx1*dy0)*dxT / (dy0*dz1 - dy1*dz0);

