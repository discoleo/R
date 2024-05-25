
# PDB Data Tools
# Exploratory Analysis


#####################

### Helper Packages & Functions

### Rpdb

# Note:
# - Function proj.line3d is available only in the version on GitHub!

# library(Rpdb)

# GitHub: 
# https://github.com/discoleo/Rpdb
#
# path = "..."
# devtools::document(path)
# devtools::load_all(path)

dist.atoms = function(p, data) {
	if(nrow(p) != 2) stop("p must be a line segment!");
	data = data$atoms[, c("x1", "x2", "x3")];
	pp = proj.line3d.matrix(data, p);
	dx = pp$x - data$x1;
	dy = pp$y - data$x2;
	dz = pp$z - data$x3;
	dd = sqrt(dx*dx + dy*dy + dz*dz);
	pp$d = dd;
	attr(pp, "d") = sqrt(sum((p[2,] - p[1,])^2));
	return(pp);
}

select.distLine = function(x, d, t, t.abs = TRUE) {
	# Note: the line-segment has t = c(0, 1);
	if(length(t) == 1) t = c(-t, t+1);
	if(t.abs) {
		d0 = attr(x, "d");
		if( ! is.null(d0)) {
			# Normalize: actual distance;
			t = t / d0;
		}
	}
	isOK = x$d <= d & (x$t >= t[1] & x$t <= t[2]);
	return(isOK);
}

residuals.pdb = function(x, filter, sep = "", ...) {
	tmp = paste(x$atoms$resname[filter], x$atoms$resid[filter], sep=sep);
	return(unique(tmp));
}

######################

### PDB Files
# setwd("...")

### PDB File: 4q02.pdb
# KRAS

x = read.pdb("4q02.pdb")

# GDP is a residue
unique(x$atoms$resname)

# Bounding Box:
isGDP = x$atoms$resname == "GDP"
bb = lapply(x$atoms[isGDP, c("x1", "x2", "x3")], range)
print(bb)

### Phosphate Groups
isP = isGDP & (x$atoms$elename %in% c("PA", "PB"))
bbP = x$atoms[isP, c("x1", "x2", "x3")];

### Ex 1:
### Distance from GDP-Channel
pp = dist.atoms(bbP, x)

# Distance from GDP-Channel
isWater = x$atoms$resname == "HOH"
isProt  = ! (isWater | isGDP)
plot(pp$d[isProt], pp$t[isProt])
points(pp$d[isWater], pp$t[isWater], col="#A032C096")
points(pp$d[isGDP], pp$t[isGDP], col="#FF0000A0")


### Neighbourhood:
isWater = x$atoms$resname == "HOH"
isClose = select.distLine(pp, d=7, t=6) & ! isGDP & ! isWater;
atoms = x$atoms[isClose, c("x1", "x2", "x3")];
resid(x, isClose)

### Plot 3D
close3d()
# Phosphate P:
points3d(bbP, col = "blue", size = 9)
points3d(x$atoms[isGDP, c("x1", "x2", "x3")], size = 6, col = "red")
points3d(atoms)

# TODO: better visualization;


### Ex 2: G-Axis
isN = isGDP & x$atoms$elename %in% c("N1", "C2", "C8")
pG  = x$atoms[isN, c("x1", "x2", "x3")];
# GDP:
close3d()
points3d(x$atoms[isGDP & ! isN, c("x1", "x2", "x3")], size = 6, col = "red")
points3d(pG, size = 6, col = "blue")
#
pG[2,] = (pG[2,] + pG[3,]) / 2;
pG = pG[c(1,2),]

pp = dist.atoms(pG, x)

### Neighbourhood:
isWater = x$atoms$resname == "HOH"
isClose = select.distLine(pp, d=7, t=7) & ! isGDP & ! isWater;
atoms = x$atoms[isClose, c("x1", "x2", "x3")];
resid(x, isClose)

### Plot 3D
# Phosphate P:
close3d()
points3d(bbP, col = "blue", size = 9)
points3d(x$atoms[isGDP, c("x1", "x2", "x3")], size = 6, col = "red")
points3d(atoms)
lines3d(pG, col = "purple", alpha = 0.75)
