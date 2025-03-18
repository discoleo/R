
### Chem 3D
# - Exploratory / Experimental code;
# - Intended to be used within Rpdb:
#   https://github.com/discoleo/Rpdb


library(rgl)

# rgl >= 1.3.1

# Math functions: are independent of rgl;
# Plot & Test functions: require rgl;


###############

### Helper

### Distances
dist.xyz = function(x, y, z) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	sqrt(diff(x)^2 + diff(y)^2 + diff(z)^2);
}
dist.l3d = function(x) {
	# TODO: process each segment of line;
	d = x[1:3] - x[4:6];
	sqrt(sum(d^2));
}
dist.p3d = function(p1, p2) {
	d = p2 - p1;
	sqrt(sum(d^2));
}

### Center
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
center.p3d = function(x) {
	dim(x) = c(3, length(x) / 3);
	apply(x, 1, mean);
}

expand.polygon3d = function(d, x, y, z, is.rel = TRUE) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	cx = mean(x); cy = mean(y); cz = mean(z);
	if( ! is.rel) {
		dx  = x - cx; dy = y - cy; dz = z - cz;
		div = sqrt(dx^2 + dy^2 + dz*2);
		d = d / div + 1; # 1 unit from centre by default;
	}
	x = d*x + (1-d)*cx;
	y = d*y + (1-d)*cy;
	z = d*z + (1-d)*cz;
	p = cbind(x=x, y=y, z=z);
	return(p);
}

### Operations

# Remove Vertex
rm.vertex.mesh = function(id, data) {
	if(length(id) > 1) {
		tmp = rm.vertex.mesh.any(id, data=data);
		return(invisible(tmp));
	}
	data$V = data$V[- id, ];
	isE = apply(data$M, 2, function(x) any(x == id));
	data$M = data$M[, ! isE];
	# Correct Edges:
	if(nrow(data$V) >= id) {
		isE = data$M > id;
		data$M[isE] = data$M[isE] - 1;
	}
	invisible(data);
}
rm.vertex.mesh.any = function(id, data) {
	data$V = data$V[- id, ];
	id  = sort(id);
	isE = apply(data$M, 2, function(x) any(x %in% id));
	data$M = data$M[, ! isE];
	for(ii in id) {
		isE = data$M >= ii;
		data$M[isE] = data$M[isE] - 1;
	}
	invisible(data);
}

### Rotations

# Rotations of Viewport
# (non-destructive rotations)

rotate.xz = function(alpha = pi/2) {
	m = matrix(c(cos(alpha),0,sin(alpha),0,
		0,1,0,0,
		-sin(alpha),0,cos(alpha),0,
		0,0,0,1), nrow = 4);
	old = par3d("userMatrix");
	par3d(userMatrix = m %*% old);
	highlevel();
	invisible(old);
}
rotate.xy = function(alpha = pi/6) {
	m = matrix(c(cos(alpha),sin(alpha),0,0,
		-sin(alpha),cos(alpha),0,0,
		0,0,1,0, 0,0,0,1), nrow = 4);
	old = par3d("userMatrix");
	par3d(userMatrix = m %*% old);
	highlevel();
	invisible(old);
}


##########################

### Orthogonal Projection/Rotation

### Projection of Point on Line
# p = point which will be projected on (x,y,z) line;
proj.line3d = function(p, x, y, z) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	dx = x[2] - x[1]; dx0 = p[1] - x[1];
	dy = y[2] - y[1]; dy0 = p[2] - y[1]; # p[2] = yP;
	dz = z[2] - z[1]; dz0 = p[3] - z[1];
	tt = (dx*dx0 + dy*dy0 + dz*dz0);
	tt = tt / (dx^2 + dy^2 + dz^2);
	t1 = 1 - tt;
	px = t1*x[1] + tt*x[2];
	py = t1*y[1] + tt*y[2];
	pz = t1*z[1] + tt*z[2];
	# TODO: compute d or not?
	# d = sqrt(sum((p - c(px, py, pz))^2));
	lst = list(P = c(px, py, pz), t = tt);
	return(lst);
}

### Projection of Point on Plane
# p = 3D Point;
# (x, y, z) = 3 points defining the plane;
proj.plane3d = function(p, x, y, z, verbose = FALSE, tol = 1E-8) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	len = c(length(x), length(y), length(z));
	stopifnot(len == 3);
	#
	id = c(1,2);
	p1 = proj.line3d(p, x[id], y[id], z[id]);
	t1 = p1$t;
	p1 = p1$P;
	pC = c(x[3], y[3], z[3]);
	pT = proj.line3d(pC, x[id], y[id], z[id]);
	tT = pT$t; dT = 1 - tT;
	if(abs(dT) < tol) {
		if(verbose) print("Case: Div 0!");
		p2 = pC * t1 + c(x[1], y[1], z[1]) * (1-t1);
	} else if(abs(t1 - 1) < tol) {
		if(verbose) print("Case: Special!");
		p2 = pC + c(x[2], y[2], z[2]) - pT$P;
	} else {
		tt = (1 - t1) / dT;
		p2 = pC * tt + c(x[2], y[2], z[2]) * (1-tt);
	}
	pp = proj.line3d(p, rbind(p1, p2));
	return(list(P = pp$P, pP12 = p1));
}


### Decompositions into Normals

eigen.lineN2 = function(p, x, y = NULL, z = NULL, normalize = TRUE) {
	N = eigen.xy3D(p, x=x, y=y, z=z, normalize=normalize);
	N = list(N1 = N$Nx, N2 = N$Ny, P = N$P);
	return(N);
}

# Both Normals: but starting in arbitrary direction;
# Out: P = projection of p0 (identical to p1);
eigen.lineAnyN2 = function(x, y = NULL, z = NULL, normalize = TRUE) {
	if(is.null(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	# Normals:
	N1 = eigen.lineAny(x=x, y=y, z=z, normalize=normalize);
	N1 = N1$N;
	p1 = c(x[1], y[1], z[1]);
	p0 = p1 + N1;
	p2 = c(x[2], y = y[2], z = z[2]);
	N2 = eigen.projPoint(p0, p1, x = p2);
	N2 = N2$N;
	#
	lst = list(N1 = N1, N2 = N2, P0 = p0, P = p1);
	return(lst);
}


### Generate 2 Normals
# (x,y,z) = coordinates of 3 points defining plane;
# p = 1st point in (x,y,z);
#     will be rotated by pi/2 around the remaining line;
### Out:
# N = normal to plane (x,y,z);
# p[out] = projected point on line;
# Note:
# - T1,T2 = rotated point (2 valid solutions);
eigen.plane = function(x, y = NULL, z = NULL, normalize = TRUE) {
	if(is.null(y)) {
		p = x[1,]; x = x[-1,];
		y = x[,2]; z = x[,3]; x = x[,1];
	} else {
		p = c(x[1], y[1], z[1]);
		x = x[-1]; y = y[-1]; z = z[-1];
	}
	# pP = Projection of P on the line;
	# pP = proj.line3d(p, x, y, z);
	# Inline variant:
	dx = x[2] - x[1]; dx0 = p[1] - x[1];
	dy = y[2] - y[1]; dy0 = p[2] - y[1];
	dz = z[2] - z[1]; dz0 = p[3] - z[1];
	tt = (dx*dx0 + dy*dy0 + dz*dz0);
	tt = tt / (dx*dx + dy*dy + dz*dz);
	t1 = 1 - tt;
	px = t1*x[1] + tt*x[2];
	py = t1*y[1] + tt*y[2];
	pz = t1*z[1] + tt*z[2];
	# Rotation:
	dx0 = p[1] - px; dx1 = x[1] - px;
	dy0 = p[2] - py; dy1 = y[1] - py;
	dz0 = p[3] - pz; dz1 = z[1] - pz;
	R = sqrt(dx0^2 + dy0^2 + dz0^2);
	dxy = dx0*dy1 - dx1*dy0;
	dxz = dx0*dz1 - dx1*dz0;
	dyz = dy0*dz1 - dy1*dz0;
	div = sqrt(dxy^2 + dxz^2 + dyz^2);
	# Note: both roots of "R" are valid;
	# TODO: check if dyz = 0, 
	dxT = dyz * R / div;
	dyT = - dxz * dxT / dyz;
	dzT =   dxy * dxT / dyz;
	# Note: - dT is also a valid solution;
	N = c(dxT, dyT, dzT);
	if(normalize) {
		# TODO: check;
		d = sqrt(dxT^2 + dyT^2 + dzT^2);
		N = N / d;
	}
	lst = list(N = N, P = c(px,py,pz));
	return(lst);
}
### Project point p on line (x,y,z)
# Returns:
# - Normals & Projection Point;
eigen.xy3D = function(p, x, y = NULL, z = NULL, normalize = TRUE) {
	if(is.null(y)) {
		xyz = x;
	} else {
		xyz = cbind(x,y,z);
	}
	xyz = rbind(p, xyz);
	N = eigen.plane(xyz, normalize=normalize);
	pP = N$P;
	Nx = pP - p;
	if(normalize) {
		Nx = Nx / sqrt(sum(Nx^2));
	}
	lst = list(Ny = N$N, Nx = Nx, P = pP);
	return(lst);
}

# Projection Point: given
# - pP = Projection point (known);
# - p0 = point which projects on pP;
# - x, y, z = forms with pP a line on which p0 is projected;
# Note: avoids recomputing pP if it is already known;
eigen.projPoint = function(p0, pP, x, y = NULL, z = NULL,
		normalize = TRUE, verbose = TRUE) {
	if(is.null(y)) {
		if(inherits(x, "matrix")) x = x[1,];
		y = x[2]; z = x[3]; x = x[1];
	}
	# Rotation:
	dx0 = p0[1] - pP[1]; dx1 = x[1] - pP[1];
	dy0 = p0[2] - pP[2]; dy1 = y[1] - pP[2];
	dz0 = p0[3] - pP[3]; dz1 = z[1] - pP[3];
	R = sqrt(dx0^2 + dy0^2 + dz0^2);
	dxy = dx0*dy1 - dx1*dy0;
	dxz = dx0*dz1 - dx1*dz0;
	dyz = dy0*dz1 - dy1*dz0;
	div = sqrt(dxy^2 + dxz^2 + dyz^2);
	# Note: both roots of "R" are valid;
	# TODO: check if dyz == 0; (already solved ?)
	NN = R / div;
	dxT =   dyz * NN;
	dyT = - dxz * NN;
	dzT =   dxy * NN;
	# Note: - dT is also a valid solution;
	N = c(dxT, dyT, dzT);
	if(normalize) {
		# TODO: may be already normalized; (check)
		d = sqrt(dxT^2 + dyT^2 + dzT^2);
		if(verbose) cat("Normalize by: ", d, "\n");
		N = N / d;
	}
	lst = list(N = N);
	return(lst);
}

# Any Normal
# x, y, z = given line;
eigen.lineAny = function(x, y = NULL, z = NULL, normalize = TRUE) {
	if(is.null(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	nx = x[2] - x[1]; ny = y[2] - y[1]; nz = z[2] - z[1];
	if(normalize) {
		d  = sqrt(nx*nx + ny*ny + nz*nz);
		nx = nx/d; ny = ny/d; nz = nz/d;
	}
	dn1 = nx - nz; dn2 = ny - nz;
	if(abs(dn2) < 1E-8) {
		# stop("Special Case: Div by 0!");
		cat("Special Case\n");
		# Normalized: if normalize == TRUE;
		# TODO: check normalization;
		fr = list(N = c(0, - ny, nz));
		return(fr);
	}
	# Arbitrary Normal:
	if(normalize) {
		div = sqrt(2 * (dn1^2 + dn2^2 - dn1*dn2));
		fx  = dn2 / div;
	} else fx = dn2;
	# fx = dn1 / div;
	fy = - fx * dn1 / dn2;
	fz = - (fx + fy);
	fr = list(N = c(fx, fy, fz));
	return(fr);
}

### Rotate by angle phi
rotate.point3d = function(p, phi, x, y = NULL, z = NULL) {
	N = eigen.xy3D(p=p, x=x, y=y, z=z, normalize = FALSE);
	pP = N$P; # Projection of p on line;
	pR = pP + cos(phi) * r * N$Nx + sin(phi) * r * N$Ny;
	return(pR);
}
### Rotate by angle phi
# N = Normals given;
rotate.point3dN2 = function(p, r, phi, N) {
	pP = N$P; # Projection of p on line;
	pR = pP + cos(phi) * r * N$Nx + sin(phi) * r * N$Ny;
	return(pR);
}

rotate.ortho3d = function(p, x, y = NULL, z = NULL) {
	if(is.null(y)) {
		xyz = x;
	} else {
		xyz = cbind(x,y,z);
	}
	xyz = rbind(p, xyz);
	N  = eigen.plane(xyz, normalize = FALSE);
	pp = N$P; px = pp[1]; py = pp[2]; pz = pp[3];
	dT = N$N; dxT = dT[1]; dyT = dT[2]; dzT = dT[3];
	# Note: - dT is also a valid solution;
	xT1 = px + dxT; yT1 = py + dyT; zT1 = pz + dzT;
	xT2 = px - dxT; yT2 = py - dyT; zT2 = pz - dzT;
	lst = list(T1 = c(xT1,yT1,zT1), T2 = c(xT2,yT2,zT2),
		P = c(px,py,pz));
	return(lst);
}


### Intersections

dist.lines3d = function(xyz, xyz0, tol = 1E-8, verbose = TRUE) {
	p1 = proj.line3d(xyz[1,], xyz0);
	d1 = sqrt(sum((xyz[1,] - p1$P)^2));
	if(d1 <= tol) {
		# Lines intersect
		if(verbose) cat("Lines intersect at P1.\n");
		return(0);
	}
	p2 = proj.line3d(xyz[2,], xyz0);
	d2 = sqrt(sum((xyz[2,] - p2$P)^2));
	if(d2 <= tol) {
		# Lines intersect
		if(verbose) cat("Lines intersect at P2.\n");
		return(0);
	}
	# Special Cases:
	d = sqrt(sum((p1$P - p2$P)^2));
	# Orthogonal lines:
	if(d <= tol) {
		# Both points project on same point;
		# Project point back:
		pO = proj.line3d(p1$P, xyz);
		dO = sqrt(sum((p1$P - pO$P)^2));
		if(dO <= tol) {
			# Lines intersect
			if(verbose) cat("Lines intersect midway.\n");
			return(0);
		}
		if(verbose) cat("Lines are orthogonal.\n");
		return(dO);
	}
	# Equidistant Projections
	if(abs(d1 - d2) <= tol) {
		pM = (p1$P + p2$P) / 2;
		pO = proj.line3d(pM, xyz);
		dO = sqrt(sum((pM - pO$P)^2));
		if(dO <= tol) {
			# Lines intersect
			if(verbose) cat("Lines intersect midway.\n");
			return(0);
		}
		# Note: lines may be parallel,
		# but is not essential for distance;
		if(verbose) cat("Special Lines but NO intersection.\n");
		return(dO);
	}
	# General Case: Minimize d
	dd = xyz[2,] - p2$P;
	dy = xyz[1,] - p1$P;
	ds = dy - dd;
	div = sum(ds*ds);
	tt  = - sum(dd * ds) / div;
	tt1 = 1 - tt;
	pL1 = tt * xyz[1,] + tt1 * xyz[2,];
	pL2 = tt * p1$P + tt1 * p2$P;
	d = sqrt(sum((pL1 - pL2)^2));
	if(d <= tol) {
		if(verbose) cat("Intersection: General case.\n");
		d = 0;
	} else if(verbose) cat("General case.\n");
	return(d);
}

# Shortest Distance between 2 lines:
# - includes points delimiting that line segment;
# - NOT to be confused with proj.line3d;
proj.lines3d = function(xyz, xyz0, tol = 1E-8) {
	p1 = proj.line3d(xyz[1,], xyz0);
	d = sqrt(sum((xyz[1,] - p$P)^2));
	if(d <= tol) {
		# Lines intersect
		lst = list(P1 = xyz[1,], p2 = p$P, d = 0,
			doIntersect = TRUE, isParallel = FALSE);
		return(lst);
	}
	stop("TODO");
	# TODO
}


################

################
### Polygons ###

# N = list with 2 Normals (N1, N2);
vertex.ngon = function(n, N, r = 1, center = c(0,0,0), phi = 0) {
	n2 = 2*n - 2;
	cs = cos(seq(0, n2, length.out = n)*pi/n + phi);
	sn = sin(seq(0, n2, length.out = n)*pi/n + phi);
	x = r*(N$N1[1]*cs + N$N2[1]*sn) + center[1];
	y = r*(N$N1[2]*cs + N$N2[2]*sn) + center[2];
	z = r*(N$N1[3]*cs + N$N2[3]*sn) + center[3];
	return(cbind(x=x, y=y, z=z));
}

############
### Ring ###

# N = list/matrix with 2 normals defining plane of ring;
# n = number of segments to partition the circle;
# length(r) = internal sub-rings;
# phi = shift of mesh;
# Note:
# - visually pleasing, but may be suboptimal
#   for engineering applications;
# - it may be better to bind rings
#   with different number of segments;
vertex.ring = function(r, N, center = c(0,0,0), n = 32, phi = 0) {
	if(length(r) == 1) {
		# TODO: optimized version
		r = c(r, 0);
	}
	#
	if(inherits(N, "list")) {
		N1 = N[[1]];
		N2 = N[[2]];
	} else if(inherits(N, "matrix")) {
		N1 = N[1,];
		N2 = N[2,];
	} else stop("Please provide the Normals!");
	pi2 = 2*pi;
	ph1 = seq(0, n-1) * pi2 / n + phi;
	ph2 = ph1 + pi/n;
	sn1 = sin(ph1); cs1 = cos(ph1);
	sn2 = sin(ph2); cs2 = cos(ph2);
	len = length(r);
	isF = rep(c(TRUE, FALSE), len %/% 2);
	if(len %% 2 == 1) isF = c(isF, TRUE);
	# Vertex
	lst = lapply(seq(len), function(id) {
		isF = isF[id];
		cs = if(isF) cs1 else cs2;
		sn = if(isF) sn1 else sn2;
		x = r[id]*(cs*N1[1] + sn*N2[1]) + center[1];
		y = r[id]*(cs*N1[2] + sn*N2[2]) + center[2];
		z = r[id]*(cs*N1[3] + sn*N2[3]) + center[3];
		cbind(x, y, z);
	});
	v = do.call(rbind, lst);
	attr(v, "dim.ring") = c(n, len);
	invisible(v);
}
mesh.ring = function(r, N, center = c(0,0,0), n = 32, phi = 0) {
	V = vertex.ring(r=r, N=N, center=center, n=n, phi=phi);
	M = mesh.ring.vertex(V);
	lst = list(V = V, M = M);
	invisible(lst);
}
mesh.ring.vertex = function(V) {
	dim = attr(V, "dim.ring");
	n = dim[1]; len = dim[2];
	# Mesh
	n2  = 2*n;
	len = len - 1;
	idM = matrix(0, nrow = 3, ncol = n2*len);
	for(id in seq(len)) {
		id0 = n2*(id - 1);
		id1 = n*(id - 1);
		id2 = id1 + n;
		idM[, id0 + seq(n)] = rbind(
			id1 + seq(n),
			id1 + c(seq(2, n), 1),
			id2 + seq(n));
		if(id %% 2 == 1) {
			idM[, id0 + seq(n)] = rbind(
				id1 + seq(n),
				id1 + c(seq(2, n), 1),
				id2 + seq(n));
			idM[, id0 + n + seq(n)] = rbind(
				id2 + seq(n),
				id2 + c(seq(2, n), 1),
				id1 + c(seq(2, n), 1));
		} else {
			idM[, id0 + seq(n)] = rbind(
				id1 + seq(n),
				id1 + c(seq(2, n), 1),
				id2 + c(seq(2, n), 1));
			idM[, id0 + n + seq(n)] = rbind(
				id2 + seq(n),
				id2 + c(seq(2, n), 1),
				id1 + seq(n));
		}
	};
	invisible(idM);
}


################
### Cylinder ###

# Note:
# - rgl::cylinder3d already implements the cylinder,
#   including proper orientation based on centre-points;
# C = Center-point from where it starts;
cylinder.line3d = function(r, x, y = NULL, z = NULL) {
	if(is.null(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	N = c(x[2] - x[1], y[2] - y[1], z[2] - z[1]);
	d = sqrt(sum(N^2));
	N = N / d;
	lst = list(R = r, L = d, N = N, C = c(x[1],y[1],z[1]));
	return(lst);
}

# Vertexes & Mesh:
mesh.cylinder = function(r, x, y = NULL, z = NULL, nL = 12, nR = 16,
		type = c("Full", "Alternating")) {
	V = mesh.vertex.cylinder(r=r, x=x, y=y, z=z, nL=nL, nR=nR);
	M = mesh.cylinder.vset(V, nL=nL, nR=nR, type=type);
	invisible(M);
}

mesh.cylinder.vset = function(V, nL = 12, nR = 16,
		type = c("Full", "Alternating")) {
	# TODO: Mesh;
	len = attr(V$V, "length");
	isOdd = len[1] > len[2];
	type = match.arg(type);
	# Note: Circle has nR+1 points;
	idS2 = len[1] * (nR + 1);
	idM = matrix(0, nrow = 3, ncol = 2*nR*len[2]);
	for(id in seq(len[2])) {
		id2 = 2*(id - 1);
		idS = id2 * nR + 1;
		idE = idS + nR - 1;
		idSlice = idS:idE;
		# Note: Circle has nR+1 points;
		idV = (id - 1)*(nR + 1) + 1;
		idV = idV:(idV + nR - 1);
		idM[, idSlice] = rbind(idV, idV + 1, idV + idS2);
		idSlice = idSlice + nR;
		idM[, idSlice] = rbind(idV + 1, idV + idS2 + 1, idV + idS2);
	}
	if(type == "Full") {
		last = len[2];
		# if(isOdd) last = last - 1;
		idM2 = matrix(0, nrow = 3, ncol = 2*nR*last);
		idS2 = idS2 - nR - 1; # previous Circle;
		for(id in seq(last)) {
			id2 = 2*(id - 1);
			idS = id2 * nR + 1;
			idE = idS + nR - 1;
			idSlice = idS:idE;
			# Note: Circle has nR+1 points;
			idV = (id - 1)*(nR + 1) + nR + 2;
			idV = idV:(idV + nR - 1);
			# print("St1"); print(idSlice)
			idM2[, idSlice] = rbind(idV + 1, idV + idS2 + 1, idV + idS2);
			idSlice = idSlice + nR;
			# print("St2"); print(idSlice)
			idM2[, idSlice] = rbind(idV, idV + 1, idV + idS2);
		}
		idM = cbind(idM, idM2);
	}
	V$M = idM;
	return(V);
}

# Vertexes:
mesh.vertex.cylinder = function(r, x, y = NULL, z = NULL, nL = 12, nR = 16) {
	if( ! is.null(y)) {
		xyz = cbind(x, y, z);
	} else xyz = x;
	# Centers
	tt = seq(0, 1, length.out = nL + 1);
	tt = cbind(1 - tt, tt);
	pC = tt %*% xyz;
	# Normals:
	N = eigen.lineAnyN2(xyz, normalize = TRUE);
	N1 = N$N1;
	N2 = N$N2;
	# Start:
	phi1 = seq(0, 2*pi, length.out = nR + 1);
	cs = cos(phi1); sn = sin(phi1);
	rx = r*(cs * N1[1] + sn * N2[1]);
	ry = r*(cs * N1[2] + sn * N2[2]);
	rz = r*(cs * N1[3] + sn * N2[3]);
	cyl1 = lapply(seq(1, nL + 1, by = 2), function(id) {
		cc = pC[id, ];
		x = rx + cc[1];
		y = ry + cc[2];
		z = rz + cc[3];
		cbind(x, y, z);
	});
	len1 = length(cyl1);
	cyl1 = do.call(rbind, cyl1);
	#
	phi2 = phi1 + pi / nR;
	cs = cos(phi2); sn = sin(phi2);
	rx = r*(cs * N1[1] + sn * N2[1]);
	ry = r*(cs * N1[2] + sn * N2[2]);
	rz = r*(cs * N1[3] + sn * N2[3]);
	cyl2 = lapply(seq(2, nL + 1, by = 2), function(id) {
		cc = pC[id, ];
		x = rx + cc[1];
		y = ry + cc[2];
		z = rz + cc[3];
		cbind(x, y, z);
	});
	len2 = length(cyl2);
	cyl2 = do.call(rbind, cyl2);
	V = rbind(cyl1, cyl2);
	attr(V, "length") = c(len1, len2);
	return(list(V = V, N1 = N1, N2 = N2));
}


# Vertexes:
# p = line segment forming the long diameter of torus;
mesh.vertex.torus = function(p, r = 1, phi = c(0, 2*pi),
		nL = 16L, nR = 16L) {
	center = (p[1,] + p[2,]) / 2;
	N1 = p[2,] - p[1,];
	R  = sqrt(sum(N1^2));
	N1 = N1 / R; R = R / 2;
	# Arbitrary N as Axis of Torus:
	N = eigen.lineAnyN2(p);
	N = list(Na = N$N1, N1 = N1, N2 = N$N2);
	v = mesh.vertex.torusN(r=r, R=R, N=N, center=center, phi=phi,
		nL=nL, nR=nR);
	invisible(v);
}
# p = central axis of torus;
# t = position of torus-centre on line p;
mesh.vertex.torusAxis = function(p, R, r = 1, t = 1/2, phi = c(0, 2*pi),
		nL = 16L, nR = 16L) {
	center = (1 - t)*p[1,] + t * p[2,];
	Na = p[2,] - p[1,];
	dd = sqrt(sum(Na^2));
	Na = Na / dd;
	# Arbitrary N as Axis of Torus:
	N = eigen.lineAnyN2(p);
	N = list(Na = Na, N1 = N$N1, N2 = N$N2);
	v = mesh.vertex.torusN(r=r, R=R, N=N, center=center, phi=phi,
		nL=nL, nR=nR);
	invisible(v);
}
mesh.vertex.torusN = function(r, R, N, center = c(0,0,0), phi = c(0, 2*pi),
		nL = 16L, nR = 16L) {
	if(inherits(N, "matrix")) N = list(Na = N[1,], N1 = N[2,], N2 = N[3,]);
	# Centers
	th = seq(phi[1], phi[2], length.out = nL);
	sc = cbind(sin(th), cos(th));
	# N2 is normalized if (N1, N2) are true normals;
	N2 = sc %*% rbind(N$N1, N$N2);
	pC = R * N2 + rep(center, each = nL);
	# Start:
	phi1 = seq(0, 2*pi, length.out = nR + 1);
	cs = cos(phi1); sn = sin(phi1);
	N1 = N$Na; rsn = r*sn;
	rx = r * cs * N1[1];
	ry = r * cs * N1[2];
	rz = r * cs * N1[3];
	cyl1 = lapply(seq(1, nL, by = 2), function(id) {
		cc = pC[id, ];
		x = rx + rsn * N2[id, 1] + cc[1];
		y = ry + rsn * N2[id, 2] + cc[2];
		z = rz + rsn * N2[id, 3] + cc[3];
		cbind(x, y, z);
	});
	len1 = length(cyl1);
	cyl1 = do.call(rbind, cyl1);
	#
	phi2 = phi1 + pi / nR;
	cs = cos(phi2); sn = sin(phi2); rsn = r*sn;
	rx = r * cs * N1[1];
	ry = r * cs * N1[2];
	rz = r * cs * N1[3];
	cyl2 = lapply(seq(2, nL, by = 2), function(id) {
		cc = pC[id, ];
		x = rx + rsn * N2[id, 1] + cc[1];
		y = ry + rsn * N2[id, 2] + cc[2];
		z = rz + rsn * N2[id, 3] + cc[3];
		cbind(x, y, z);
	});
	len2 = length(cyl2);
	cyl2 = do.call(rbind, cyl2);
	V = rbind(cyl1, cyl2);
	attr(V, "length") = c(len1, len2);
	return(list(V = V, Na = N$Na, N1 = N$N1, N2 = N$N2));
}


### Cylinder: Diagonal Section
# - very basic implementation;

# p = 2 points defining central axis;
cylinder.section.p2 = function(px, p, r,
		type = c("Alternating", "Simple", "Ellipse"),
		nC = 32, nL = 17, phi = 0) {
	N = eigen.lineN2(px, p);
	lst = cylinder.section(p=p, r=r, N=N, type=type,
		nC=nC, nL=nL, phi=phi);
	invisible(lst);
}

# p = 2 points defining central axis;
cylinder.section = function(p, r, N,
		type = c("Alternating", "Simple", "Ellipse"),
		nC = 32, nL = 17, phi = 0) {
	type = match.arg(type);
	if(type == "Simple") {
		lst = cylinder.section.simple(p, r=r, N=N, nC=nC, nL=nL, phi=phi);
		return(lst);
	} else if(type == "Ellipse") {
		lst = cylinder.section.ellipse(p, r=r, N=N, nC=nC, nL=nL, phi=phi);
		return(lst);
	}
	# Alternating
	lst = cylinder.section.alternating(p, r=r, N=N, nC=nC, nL=nL, phi=phi);
	return(lst);
}
cylinder.section.ellipse = function(p, r, N, nC = 32, nL = 17, phi = 0) {
	# Half of Cylinder = Centre of Diagonal;
	p[2,] = (p[1,] + p[2,])/2;
	t0 = seq(0, 1, length.out = nL);
	ct = cbind(1 - t0, t0) %*% p;
	dp = p[2,] - p[1,];
	d2 = sum(dp^2);
	rt = sqrt(r^2 + d2*t0^2);
	Nt = lapply(seq_along(t0), function(id) {
		(t0[id] * dp + r*N$N1) / rt[id];
	});
	Nt = do.call(rbind, Nt);
	tc = seq(0, 2*pi, length.out = nC) + phi;
	pp = lapply(seq_along(t0), function(id) {
		x = rt[id] * cos(tc) * Nt[id,1] + r * sin(tc) * N$N2[1] + ct[id,1];
		y = rt[id] * cos(tc) * Nt[id,2] + r * sin(tc) * N$N2[2] + ct[id,2];
		z = rt[id] * cos(tc) * Nt[id,3] + r * sin(tc) * N$N2[3] + ct[id,3];
		cbind(x,y,z);
	});
	pp = do.call(rbind, pp);
	invisible(pp);
}
cylinder.section.simple = function(p, r, N,
		nC = 32, nL = 17, phi = 0) {
	dp = p[2,] - p[1,];
	d2 = sum(dp^2);
	rt = sqrt(r^2 + d2/4);
	# Nt = (dp/2 + r*N$N1) / rt;
	# Intersection: Circles w Ellipse
	t0 = seq(0, 1, length.out = nL);
	yr = 2 * sqrt(t0*(1 - t0));
	th = acos(yr);
	th[t0 <= 1/2] = th[t0 <= 1/2] + pi/2;
	th[t0 > 1/2]  = pi/2 - th[t0 > 1/2];
	ct = cbind(1 - t0, t0) %*% p;
	# Note: tc can be made with alternating phase;
	tc = seq(0, 2*pi, length.out = nC);
	pp = lapply(seq_along(t0), function(id) {
		tc = tc[tc < th[id]];
		# add boundary point:
		tc = c(tc, th[id]);
		tc = c(tc, - tc) + phi;
		x = r * (cos(tc) * N$N1[1] + sin(tc) * N$N2[1]) + ct[id,1];
		y = r * (cos(tc) * N$N1[2] + sin(tc) * N$N2[2]) + ct[id,2];
		z = r * (cos(tc) * N$N1[3] + sin(tc) * N$N2[3]) + ct[id,3];
		cbind(x,y,z);
	});
	pp = do.call(rbind, pp);
	invisible(pp);
}
cylinder.section.alternating = function(p, r, N,
		nC = 32, nL = 17, phi = 0) {
	dp = p[2,] - p[1,];
	d2 = sum(dp^2);
	rt = sqrt(r^2 + d2/4);
	# Nt = (dp/2 + r*N$N1) / rt;
	# Intersection: Circles w Ellipse
	t0 = seq(0, 1, length.out = nL);
	yr = 2 * sqrt(t0*(1 - t0));
	th = acos(yr);
	th[t0 <= 1/2] = th[t0 <= 1/2] + pi/2;
	th[t0 > 1/2]  = pi/2 - th[t0 > 1/2];
	ct = cbind(1 - t0, t0) %*% p;
	# Alternating phase;
	tc  = seq(0, 2*pi, length.out = nC);
	pp1 = lapply(seq(1, nL, by=2), function(id) {
		tc = tc[tc < th[id]];
		# add boundary point:
		tc = c(tc, th[id]);
		tc = c(tc, - tc) + phi;
		x = r * (cos(tc) * N$N1[1] + sin(tc) * N$N2[1]) + ct[id,1];
		y = r * (cos(tc) * N$N1[2] + sin(tc) * N$N2[2]) + ct[id,2];
		z = r * (cos(tc) * N$N1[3] + sin(tc) * N$N2[3]) + ct[id,3];
		cbind(x,y,z);
	});
	pp1 = do.call(rbind, pp1);
	# Phased:
	tc  = tc + pi / nC;
	pp2 = lapply(seq(2, nL, by=2), function(id) {
		tc = tc[tc < th[id]];
		# add boundary point:
		tc = c(tc, th[id]);
		tc = c(tc, - tc) + phi;
		x = r * (cos(tc) * N$N1[1] + sin(tc) * N$N2[1]) + ct[id,1];
		y = r * (cos(tc) * N$N1[2] + sin(tc) * N$N2[2]) + ct[id,2];
		z = r * (cos(tc) * N$N1[3] + sin(tc) * N$N2[3]) + ct[id,3];
		cbind(x,y,z);
	});
	pp2 = do.call(rbind, pp2);
	pp  = rbind(pp1, pp2);
	invisible(pp);
}


############
### Cone ###

cone.vertex.alternating = function(r, p, nL = 16, nC = 32, phi = 0) {
	# Arbitrary N to Axis of Cone:
	N = eigen.lineAnyN2(p);
	# Axis:
	Na = p[2,] - p[1,];
	dd = sqrt(sum(Na^2));
	Na = Na / dd;
	# Centres:
	tp = seq(0, 1, length.out = nL);
	ti = 1 - tp;
	ct = cbind(tp, ti) %*% p;
	# Circles
	tc = seq(0, 2*pi, length.out = nC) + phi;
	sc = cbind(cos(tc), sin(tc));
	v1 = lapply(seq(1, nL, by=2), function(id) {
		r = r * ti[id];
		x = r * (sc[,1] * N$N1[1] + sc[,2] * N$N2[1]) + ct[id,1];
		y = r * (sc[,1] * N$N1[2] + sc[,2] * N$N2[2]) + ct[id,2];
		z = r * (sc[,1] * N$N1[3] + sc[,2] * N$N2[3]) + ct[id,3];
		cbind(x,y,z);
	});
	v1 = do.call(rbind, v1);
	#
	tc = tc + pi / nC;
	sc = cbind(cos(tc), sin(tc));
	v2 = lapply(seq(2, nL, by=2), function(id) {
		r = r * ti[id];
		x = r * (sc[,1] * N$N1[1] + sc[,2] * N$N2[1]) + ct[id,1];
		y = r * (sc[,1] * N$N1[2] + sc[,2] * N$N2[2]) + ct[id,2];
		z = r * (sc[,1] * N$N1[3] + sc[,2] * N$N2[3]) + ct[id,3];
		cbind(x,y,z);
	});
	v2 = do.call(rbind, v2);
	vv = rbind(v1, v2);
	invisible(vv);
}

cone.vertex.adaptive = function(r, p, nL = 16, nC = 32, phi = 0,
		options = list(min = 6, limit = 6, shift.fr = 1)) {
	# Arbitrary N to Axis of Cone:
	N = eigen.lineAnyN2(p);
	# Axis:
	Na = p[2,] - p[1,];
	dd = sqrt(sum(Na^2));
	Na = Na / dd;
	# Centres:
	tp = seq(0, 1, length.out = nL);
	ti = 1 - tp;
	ct = cbind(tp, ti) %*% p;
	# Circles
	nt = ceiling(nC * ti);
	# Very small circles:
	nt[nt < options$limit] = options$min;
	nt[length(nt)] = 1;
	#
	dShift = options$shift.fr;
	vv = lapply(seq_along(tp), function(id) {
		tc = seq(0, 2*pi, length.out = nt[id]) + phi;
		if(id %% 2 == 0) {
			# TODO: best phase shift ???
			tc = tc + pi / (nt[id] + dShift);
		}
		sc = cbind(cos(tc), sin(tc));
		r = r * ti[id];
		x = r * (sc[,1] * N$N1[1] + sc[,2] * N$N2[1]) + ct[id,1];
		y = r * (sc[,1] * N$N1[2] + sc[,2] * N$N2[2]) + ct[id,2];
		z = r * (sc[,1] * N$N1[3] + sc[,2] * N$N2[3]) + ct[id,3];
		cbind(x,y,z);
	});
	vv = do.call(rbind, vv);
	invisible(vv);
}

cone.vertex.simple = function(r, p, nL = 16, nC = 32, phi = 0) {
	# Arbitrary N to Axis of Cone:
	N = eigen.lineAnyN2(p);
	# Axis:
	Na = p[2,] - p[1,];
	dd = sqrt(sum(Na^2));
	Na = Na / dd;
	# Centres:
	tp = seq(0, 1, length.out = nL);
	ti = 1 - tp;
	ct = cbind(tp, ti) %*% p;
	# Circles
	tc = seq(0, 2*pi, length.out = nC) + phi;
	sc = cbind(cos(tc), sin(tc));
	vv = lapply(seq_along(tp), function(id) {
		r = r * ti[id];
		x = r * (sc[,1] * N$N1[1] + sc[,2] * N$N2[1]) + ct[id,1];
		y = r * (sc[,1] * N$N1[2] + sc[,2] * N$N2[2]) + ct[id,2];
		z = r * (sc[,1] * N$N1[3] + sc[,2] * N$N2[3]) + ct[id,3];
		cbind(x,y,z);
	});
	vv = do.call(rbind, vv);
	invisible(vv);
}


#############
### Helix ###

helix = function(r, x, y = NULL, z = NULL, d = 0.25) {
	if( ! is.null(y)) {
		xyz = cbind(x, y, z);
	} else xyz = x;
	# Normals:
	N = eigen.lineAnyN2(xyz, normalize = TRUE);
	N1 = N$N1;
	N2 = N$N2;
	N0 = xyz[2,] - xyz[1,];
	d  = sqrt(sum(N0^2));
	N0 = N0 / d;
	# TODO: better concept;
	Ns = N2 + N0; # Helix Step
	p3 = rbind(xyz[1,] + Ns, xyz[1,], xyz[1,] + N1 + Ns);
	Nw = eigen.plane(p3);
	N$Nw = Nw$N; N$N0 = N0; N$d = d;
	return(N);
}

# Simple Ribbon:
helix.ribbon = function(p, N, helix = 4, r = 1, col = NULL, alpha = 1,
		type = c("Simple", "Screw"), dashed = TRUE,
		n.points = 32, dth = pi/7, scale.screw = 1.5) {
	np = helix * n.points;
	type = match.arg(type);
	hasC = FALSE;
	if(is.null(dim(p))) {
		p1 = p;
	} else {
		p1 = p[1,];
		if(dim(p)[1] > 1) hasC = TRUE;
	}
	phi0 = seq(0, np) * 2*pi / n.points;
	dth  = c(0, -dth, dth);
	if(type == "Screw") {
		smR = r * cos(dth[3] / scale.screw);
		r = c(r, smR, smR);
	} else {
		r = c(r,r,r);
	}
	if(dashed) {
		# Note: more efficient than drawing each line-segment individually;
		alpha = rep(c(alpha, 0), np %/% 2);
		if(np %% 2 == 1) alpha = c(alpha, alpha[1]);
	}
	pp = lapply(seq(3), function(idRibbons) {
		th = dth[idRibbons];
		rr = r[idRibbons];
		pp = lapply(seq(np + 1), function(id) {
			phi = phi0[id] + th;
			pp = rr * (cos(phi)*N$N2 - sin(phi)*N$N1) +
				+ N$N0*N$d*(id - 1)/np + p1;
			return(pp);
		})
		pp = do.call(rbind, pp);
		list(V = pp, col=col, alpha=alpha);
	});
	pp = list(Ribbon = pp);
	if(hasC) pp$P = p;
	invisible(pp);
}

plot.helix.ribbon = function(v, col = NULL) {
	isColNull = is.null(col);
	if( ! is.null(v$P)) {
		lines3d(v$P);
	}
	vR = v$Ribbon;
	lapply(vR, function(v) {
		colRib = v$col;
		if(is.null(colRib)) {
			if(isColNull) {
				lines3d(v$V, alpha = v$alpha);
			} else {
				lines3d(v$V, alpha = v$alpha, col=col);
			}
		} else {
			lines3d(v$V, alpha = v$alpha, col = colRib);
		}
	});
	invisible();
}


###############

### Tetrahedron
Th4 = function(r = 1, center = c(0,0,0), phi = 0) {
	rT = r; r = r * 2 * sqrt(2) / 3;
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
	return(list(V = pV, C = pC, R = rT));
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
		id = if(id.base == 1) c(1,2,3)
			else if(id.base == 2) c(1,2,4)
			else if(id.base == 2) c(1,3,4) else c(2,3,4);
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

###
# xyz    = c(a,b,c)
# center = center of tetrahedron;
Th4 = function(xyz, center = c(0,0,0), r = 1) {
	# Eq: a*(x-x0) + b*(y-y0) + c*(z-z0) = 0
	# TODO
}


###################

###################
### Icosahedron ###

mesh.icosa.vset = function(v) {
	m = cbind(
		rbind(1:5,  c(2:5, 1),   6:10),
		rbind(1:5,    6:10,   c(10, 6:9)),
		rbind(1:5,  c(2:5, 1),  11),
		rbind(6:10, c(7:10, 6), 12)
	);
	v = list(V = v, M = m);
	invisible(v);
}


# p  = points defining Axis;
# d  = side-length; default = based on ||p||;
# dH = distance between 2 pentagon-planes;
# dV = distance between vertex and pentagon;
# Note: changing dH or dV generates a non-regular Icosahedron;
icosa = function(p, d = NULL, dH = NULL, dV = NULL, phi = 0, detailed = FALSE) {
	### Formulas:
	# Rp = r of poly(5);
	# d  = 2*Rp*sin(pi/5);
	# d0 = dist(topV, bottomV);
	# dd = dist(Vertex, Diagonal);
	d0  = dist.xyz(p);
	thi = 1/tan(pi/5);
	dd  = (1/sin(pi/5) - thi) / 2;
	dH1 = sqrt(3/4 - dd^2);
	dV1 = sqrt(3/4 - thi^2/4);
	isNH = is.null(dH);
	isNV = is.null(dV);
	# TODO: check concept & combinations;
	# - may require different concept if some of the values are given;
	if(isNH) dH = dH1;
	if(isNV) dV = dV1;
	if(is.null(d)) {
		d = d0 / (dH + 2*dV);
	}
	if(isNH) dH = dH * d;
	if(isNV) dV = dV * d;
	# Centres:
	t1 = dV / d0;
	t2 = 1 - t1;
	c1 = t2 * p[1,] + t1 * p[2,];
	c2 = t1 * p[1,] + t2 * p[2,];
	# Pentagons:
	N  = eigen.lineAnyN2(p);
	rp = d / (2*sin(pi/5));
	v1 = vertex.ngon(5, N=N, r = rp, center = c1, phi = phi);
	v2 = vertex.ngon(5, N=N, r = rp, center = c2, phi = phi + pi/5);
	v  = rbind(v1, v2, p);
	if(detailed) {
		attr(v, "details") = list(t = c(t1, t2), C = rbind(c1, c2));
	}
	return(v);
}

# r = radius of circumscribed sphere;
# N = Normal of main (user defined) axis;
vertex.icosa.r = function(r, N, center = c(0,0,0), phi = 0, detailed = TRUE) {
	lst = math.icosa.r(r);
	rp  = lst$dS / (2*sin(pi/5));
	# Centres:
	p1 = center - r*N;
	p2 = center + r*N;
	t1 = lst$t1; t2 = lst$t2;
	c1 = t2 * p1 + t1 * p2;
	c2 = t1 * p1 + t2 * p2;
	# Pentagons:
	N2 = eigen.lineAnyN2(rbind(c(0,0,0), N));
	v1 = vertex.ngon(5, N = N2, r = rp, center = c1, phi = phi);
	v2 = vertex.ngon(5, N = N2, r = rp, center = c2, phi = phi + pi/5);
	v  = rbind(v1, v2, p1, p2);
	if(detailed) {
		attr(v, "details") = list(t = c(t1, t2), C = rbind(c1, c2));
	}
	return(v);
}
math.icosa.r = function(r) {
	thi = 1/tan(pi/5);
	dd =  1/sin(pi/5) - thi;
	dH = sqrt(3 - dd^2) / 2;
	dV = sqrt(3 - thi^2) / 2;
	r2 = 2*r;
	d  = r2 / (dH + 2*dV);
	dH = dH * d;
	dV = dV * d;
	# Centres:
	t1 = dV / r2;
	t2 = 1 - t1;
	# Result:
	lst = list(R = r, dS = d, dH=dH, dV=dV, t1=t1, t2=t2);
}

### Plot:
plot.bb.icosa = function(v, p = NULL,
		col.line = "blue", col.p = c("red", "orange"), size = 5) {
	# Main (custom) Axis:
	if(is.null(p)) p = v[c(11,12), ];
	points3d(v);
	lines3d(v[c(1:5, 1),]);
	lines3d(v[c(6:10, 6),]);
	lines3d(p, col = col.line);
	points3d(p, col = col.p[[1]], size=size);
	C = attr(v, "details")$C;
	if( ! is.null(C)) {
		points3d(C, col = col.p[[2]], size=size);
	}
	invisible(v);
}


# TODO: regular anti-prism;


##################

### Torus: Cyclo-6
# Poly-cyclic compound around a cylinder;
cyclo.cylinder = function(n, R = 3) {
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
	obj = lapply(seq(n), function(id) {
		cbind(
			x = asV(x1, xB, xS, id=id),
			y = asV(y1, yB, yS, id=id),
			z = c(0, dz, dz + ds, dz, 0, - ds, 0));
	});
	invisible(obj);
}


###################
###################

### Spheres Inside a Circle
# - 3D version: inside larger circle of given radius;
# - Initial code in pseudo-package BioShapes, see:
#   Graphic.Helper.Circles.Chains.R;
circles.InFixedCircle = function(n, r, N, center = c(0,0,0), phi=0) {
	if(inherits(N, "list")) {
		N1 = N[[1]]; N2 = N[[2]];
	} else if(inherits(N, "matrix")) {
		N1 = N[1,]; N2 = N[2,];
	} else stop("Please provide valid Normals!");
	r1 = r / (1 + 1/sin(pi/n));
	R  = r - r1;
	th = seq(0, n-1) * 2*pi / n + phi;
	sn = sin(th); cs = cos(th);
	x = R*(cs*N1[1] + sn*N2[1]) + center[1];
	y = R*(cs*N1[2] + sn*N2[2]) + center[2];
	z = R*(cs*N1[3] + sn*N2[3]) + center[3];
	xyz = list(center = cbind(x,y,z), r = r1, R = R);
	return(xyz);
}

### Outer Bearing
bearing = function(r, N, center = c(0,0,0), n = 15,
		h.scale = 0.75, phi = 0) {
	N3 = N[3,];
	bb = circles.InFixedCircle(n=n, r=r, N=N, center=center, phi=phi);
	# Inner wall:
	rb = bb$r;
	ri = bb$R - rb; # inner radius
	pp = rbind(center - rb*N3, center + rb*N3);
	bt = mesh.cylinder(r = ri, pp);
	# Lateral walls:
	hh = (1+h.scale) * rb;
	rw = if(h.scale > 0) c(ri + hh, ri + rb, ri) else c(ri + hh, ri);
	w1 = mesh.ring(r = rw, N, center = pp[1,]);
	w2 = mesh.ring(r = rw, N, center = pp[2,]);
	# Outer wall:
	# TODO
	# Bearing:
	lst = list(Balls = bb, Base = bt, Wall1 = w1, Wall2 = w2);
	invisible(lst);
}

### Plot:
plot.spheres.chain = function(x, ...) {
	cc  = x$center;
	len = nrow(cc);
	for(id in seq(len)) {
		spheres3d(cc[id,], y = NULL, z = NULL, radius = x$r, ...);
	}
	invisible();
}

plot.bearing = function(x, col = c("#323232", "#969696", "#A0A0A0"),
		alpha = c(1, 0.8, 0.6, 0.6), ...) {
	if(length(col) == 3) col = c(col, col[[3]]);
	plot.spheres.chain(x$Balls, col = col[[1]], alpha = alpha[[1]]);
	mCyl = mesh3d(vertices = t(x$Base$V), triangles = x$Base$M)
	shade3d(mCyl, col = col[[2]], alpha = alpha[[2]]);
	# Lateral Walls:
	mW1 = mesh3d(vertices = t(x$Wall1$V), triangles = x$Wall1$M);
	mW2 = mesh3d(vertices = t(x$Wall2$V), triangles = x$Wall2$M);
	shade3d(mW1, col = col[[3]], alpha = alpha[[3]]);
	shade3d(mW2, col = col[[4]], alpha = alpha[[4]]);
}

