

library(rgl)

# rgl >= 1.3.1


###############

### Helper

dist.xyz = function(x, y, z) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	sqrt(diff(x)^2 + diff(y)^2 + diff(z)^2);
}

center.xyz = function(x, y, z) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	len = length(x);
	x = sum(x) / len;
	y = sum(y) / len;
	z = sum(z) / len;
	return(c(x,y,z));
}

###############

### Tetrahedron
Th4 = function(r = 1, center = c(0,0,0), phi = 0) {
	# Base Triangle:
	pB = if(phi == 0) {
			r2 = - r/2; r3 = - r2 * sqrt(3);
			cbind(c(r,r2,r2) + center[1], c(0,r3,-r3) + center[2], center[3]);
		} else {
			phi3 = c(0, 2*pi/3, 4*pi/3) + phi;
			x = r * cos(phi3) + center[1];
			y = r * sin(phi3) + center[2];
			cbind(x, y, center[3]);
		}
	pT = c(0, 0, r * sqrt(2)) + center; # Top
	pV = rbind(pB, pT);
	pC = c(0, 0, r/4 * sqrt(2)) + center; # Center
	return(list(V = pV, C = pC, R = r));
}

plot.Th4 = function(v, lwd = 2, id.base = NULL,
		fill.plane = "#A0A0A064", N = 20, alpha = 0.2) {
	# Lines from Center to Vertices:
	pV = v$V;
	pC = v$C;
	for(i in seq(4)) {
		tmp = pV[i,];
		lines3d(c(tmp[1], pC[1]), c(tmp[2], pC[2]), c(tmp[3], pC[3]), lwd=lwd);
	}
	# Draw Base-plane:
	if( ! is.null(id.base)) {
		# TODO:
		id = c(1,2,3);
		pB = as.data.frame(pV[id, ]);
		names(pB) = c("x","y","z");
		### Transparency: set with alpha!
		polygon3d(pB, col = fill.plane, fill = TRUE, alpha = alpha);
		# fake transparency:
		trf = function(x, tt) {
			x1 = tt * x[1] + (1 - tt) * x[2];
			x2 = tt * x[1] + (1 - tt) * x[3];
			x  = c(x1, x2);
		}
		# Hashlines in Base-Plane:
		if(N > 0)
		for(tt in seq(0, N - 1)/N) {
			x = trf(pB$x, tt); y = trf(pB$y, tt);
			z = trf(pB$z, tt);
			lines3d(x, y, z, lwd = lwd, col = fill.plane);
		}
	}
}

### Basic Tetrahedron
# Note: transparency works only by specifying alpha!
Th4.base = function(lwd = 2, fill.plane = "#A0A0A064", N = 20, alpha = 0.2) {
	# Base Triangle:
	pB = list(c(1,0,0), c(-1/2, sqrt(3)/2, 0), c(-1/2, -sqrt(3)/2, 0))
	pT = c(0, 0, sqrt(2)); # Top
	pC = c(0, 0, 1/4 * sqrt(2)); # Center
	# Lines from Center to Vertices:
	for(i in seq(3)) {
		tmp = pB[[i]];
		lines3d(c(tmp[1], pC[1]), c(tmp[2], pC[2]), c(tmp[3], pC[3]), lwd=lwd);
	}
	tmp = pT;
	lines3d(c(tmp[1], pC[1]), c(tmp[2], pC[2]), c(tmp[3], pC[3]), lwd=lwd);
	# Draw Base-plane:
	if( ! is.null(fill.plane)) {
		pB = data.frame(do.call("rbind", pB));
		names(pB) = c("x","y","z");
		# Bug in old rgl:::triangulateSimple
		# polygon3d(pB, col = fill.plane);
		# polygon3d(rbind(pB, pB[1,]), ...);
		### Transparency: set with alpha!
		polygon3d(pB, col = fill.plane, fill = TRUE, alpha = alpha);
		# fake transparency:
		trf = function(x, tt) {
			x1 = tt * x[1] + (1 - tt) * x[2];
			x2 = tt * x[1] + (1 - tt) * x[3];
			x  = c(x1, x2);
		}
		# Hashlines in Base-Plane:
		if(N > 0)
		for(tt in seq(0, N - 1)/N) {
			x = trf(pB$x, tt); y = trf(pB$y, tt); z = c(0, 0);
			lines3d(x, y, z, lwd = lwd, col = fill.plane);
		}
	}
}

close3d()
open3d()
Th4.base(N = 32)


###
# xyz    = c(a,b,c)
# center = center of tetrahedron;
Th4 = function(xyz, center = c(0,0,0), r = 1) {
	# Eq: a*(x-x0) + b*(y-y0) + c*(z-z0) = 0
	# TODO
}


##################

### Torus: Cyclo-6
cycloCylinder = function(n, R = 3) {
	# Th:
	# th = pi - atan(4); cs = - 1/sqrt(17);
	# hb = dh / (1 + 1/sqrt(17));
	dz  = R * sin(pi/n) / (1 + 1/sqrt(17));
	# TODO:
	ds  = dz / 1.2;
	x1 = R * cos(2*pi*seq(n)/n); x1 = c(x1, x1[1]);
	y1 = R * sin(2*pi*seq(n)/n); y1 = c(y1, y1[1]);
	xB = 1.25 * R * cos(2*pi*seq(n)/n + pi/n);
	yB = 1.25 * R * sin(2*pi*seq(n)/n + pi/n);
	xS = 0.75 * R * cos(2*pi*seq(n)/n + pi/n);
	yS = 0.75 * R * sin(2*pi*seq(n)/n + pi/n);
	asV = function(x, xB, xS, id) {
		c(x[id], x[id], xB[id], x[id + 1], x[id + 1], xS[id], x[id]);
	}
l	tmp = lapply(seq(n), function(id) {
		lines3d(
			x = asV(x1, xB, xS, id=id),
			y = asV(y1, yB, yS, id=id),
			z = c(0, dz, dz + ds, dz, 0, - ds));
	})
	invisible()
}


# Top NMR specialists have gathered in an emergency
# to clarify which conformation is adopted: chair or boat?
close3d()
open3d()
cycloCylinder(19)

