

### Oncology: Mathematical Models


# library(BioShapes)
# see https://github.com/discoleo/BioShapes


### Save Images
# setwd(...)


##################

### Tumor Invasion

### Gatenby–Gawlinski Model:

### Ref:
# 1. Gatenby RA, Gawlinski ET. A reaction-diffusion model of cancer invasion.
#    Cancer Res. (1996); 56(24):5745–5753;
#
### Extension of the initial Gatenby-Gawlinski model:
# 2. McGillen JB, Gaffney EA, Martin NK, Maini PK.
#    A general reaction–diffusion model of acidity in cancer invasion.
#    J. Math. Biol. (2014); 68:1199–1224;
#    http://doi.org/10.1007/s00285-013-0665-7

### Time: t = 0
plot.diagram.t0 = function(dL = 0.0125, asp = 4/3, new = TRUE,
		cex = 1.6, cex.axis = cex - 0.1) {
	if(new) plot.base(xlim = c(-1, 1), ylim = c(0, 1.5),
		cex.axis = cex.axis, asp = asp);
	
	lines(c(-1,0), c(1,0), lwd=3, col = "red")
	lines(c(-1,0,1), c(0,1,1), lwd=2, col = "#323232")
	abline(v = 0, lty=2, col = "green")
	
	polygon(c(0,1,1,0), c(0,0,1,1), density = 8, col = "#64646480")
	text(0.5, 0.5, "Normal\ntissue", cex=2)
	# Tumor:
	arrowSimple(c(-0.5, -0.8), c(1.3, 0.9), col = "red",
		d = -0.1, d.lines = c(-dL,dL), scale = 1/asp)
	text(-0.5, 1.3, "Density of Tumor Cells", adj = c(-0.05, -0.05),
		cex = cex, col = "red")
	# Boundary:
	arrowSimple(c(0.3, 0.025), c(1.2, 1.2), col = "green",
		d = -0.1, d.lines = c(-dL,dL), scale = 1/asp)
	text(0.35, 1.2, "Boundary:", adj = c(0, 0.5),
		cex = cex, col = "green")
	text(0, 1.1, "Tumor - Normal Tissue", adj = c(-0.05, 0.5),
		cex = cex, col = "green")
	# Title:
	text(-0.8, 1.45, "Invasion Model:    t = 0", adj = c(0, 0.5),
		cex = cex + 0.2, col = "blue")
}

# png(file = "Model.Gatenby.T0.png")

plot.diagram.t0()

# dev.off()


###################

###################
### Sub-Populations


### Ref:
# 2. Martin NK, Gaffney EA, Gatenby RA, Maini PK.
#    Tumour-stromal interactions in acid-mediated invasion: a mathematical model.
#    J Theor Biol. 2010 Dec 7;267(3):461-70.
#    http://doi.org/10.1016/j.jtbi.2010.08.028.
#    Epub 2010 Sep 8. PMID: 20816684; PMCID: PMC3005191.
#  - model with ECM-degradation;
# 3. Strobl MAR, Krause AL, Damaghi M, Gillies R, Anderson ARA, Maini PK.
#    Mix and Match: Phenotypic Coexistence as a Key Facilitator of Cancer Invasion.
#    Bull Math Biol. 2020 Jan 17;82(1):15.
#    http://doi.org/10.1007/s11538-019-00675-0
#    PMID: 31953602; PMCID: PMC6968991.
#  - model with 2 sub-populations of tumor cells;

# TODO:

# png("img.Tumor.States.png")

lim = c(-1, 10)
plot.base(xlim = lim, ylim = lim, axt = NULL)
# plot.base(xlim = lim, ylim = lim)
text(
	c(1, 7, 1, 7, 4),
	c(8, 8, 5, 5, 2),
	c("Acid", "Stroma", "Ta", "Tm", "ECM"), cex=2)
lwd = 2;
arrowT(c(2,5.5), c(8,8), lwd=lwd)
arrowT(c(6,2), c(7.5, 5.5), lwd=lwd)
arrowT(c(7,7), c(7.5, 5.5), lwd=lwd)
arrowTriangle(c(1,1), c(5.5, 7.5), lwd=lwd, col = "#64B032", fill = "#64B032")
# Ta |---| Tm
arrowT(c(2,6), c(5,5) + 0.2, lwd=lwd)
arrowT(c(6,2), c(5,5) - 0.2, lwd=lwd)
# ECM
arrowT(c(3.5, 1.5), c(2.5, 4.5), lwd=lwd)
arrowT(c(4.5, 6.5), c(2.5, 4.5), lwd=lwd)
arrowT(c(7.0, 5.0), c(4.5, 2.5), lwd=lwd)

# dev.off()

