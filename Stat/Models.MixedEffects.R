#######################
##
## Statistical Series:
## Mixed Effects Models
##
## Leonard Mada

### Mixed Effects Models
# - Stats Lecture for Medical Residents:
#   Plots & Images;


####################

library(lattice)
library(nlme)


### Helper Functions

# Note:
# - it proved easier to modify the generated Trellis-object,
#   then set the text & symbol size;

# Images for Presentations:
# - Text: larger size required;
setCEX = function(x, cex = 1.2, text.size = 16, points.size = 12,
		cex.points = 1, col.line = NULL) {
	hasLegend = ! is.null(x$legend);
	if(hasLegend) x$legend$top$args$key$text$cex = cex;
	# Axis:
	x$x.scales$cex = c(cex, cex);
	x$y.scales$cex = c(cex, cex);
	# Panels:
	# $regions$col; $superpose.line$col; $plot.symbol$col
	# x$panel.args.common
	# x$panel.args.common$col.line
	# x$panel.args.common$cex = cex.points;
	if(is.null(x$par.settings)) x$par.settings = list();
	fsz = list(text = text.size, points = points.size);
	x$par.settings$fontsize = fsz;
	if(! is.null(col.line)) {
		# if(is.null(x$par.settings$superpose.line)) x$par.settings$superpose.line = list();
		# x$par.settings$superpose.line$lwd = 4;
		# x$par.settings$superpose.line$col = rep(col.line, 7);
		# x$par.settings$plot.symbol = list(col = col.line)
		x$par.settings$plot.line = list(col = col.line, lwd = 1.75);
	}
	return(x);
}

######################
######################

# setwd(...)

IMG = FALSE;

# Note:
# - Sizes are intended for images included
#   in the presentation; (not for real-time display)

######################

### Example: Orthodont

# Distance: Pituitary - Pterygo-maxillary fissure

###
xAll = Orthodont;
attr(xAll, "labels")$y = "Distance\nPituitary - Pterygomax. fissure";
xF = xAll[xAll$Sex == "Female", ];


### Plot Dataset
if(IMG) png(file = "img/LM.Mixed.Orthodont.png", width = 1200, height = 760)

tmp = plot(xAll, innerGroups = ~ Sex)
setCEX(tmp, text.size = 20)

dev.off()


### Subset: Females
if(IMG) png(file = "img/LM.Mixed.Orthodont.F.png", width = 800, height = 380)

tmp = plot(xF, col.symbol = 2)
setCEX(tmp, col.line = "#F0A232")

dev.off()


### Stats:
fm1 = lme(distance ~ age, data = Orthodont) # random is ~ age
fm2 = lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
summary(fm1)
summary(fm2)

# TODO

#######################

### Example: BodyWeight
# Body Weight of Rats

colDiet = c(2:4)


head(BodyWeight)
table(BodyWeight$Diet)


plot(BodyWeight, innerGroups = ~ Diet, col = 2:4,
	panel = function(x, y, groups, subscripts, ...) {
		panel.xyplot(x, y, ..., groups = groups, subscripts = subscripts, )
	})


plot(BodyWeight, innerGroups = ~Diet)


###
rat.fit = lme(weight ~ Time + Diet, data = BodyWeight, random = ~ Time|Rat)

# Note: wastes panel-space;
plot(rat.fit, weight ~ Time | Rat + Diet, groups = ~Diet, col = 2:4,
     panel = function(x, y, groups, subscripts, ...) {
  panel.xyplot(x, y, ..., groups = groups, subscripts = subscripts, )
})

plot(rat.fit, weight ~ Time | Rat + Diet)

plot(rat.fit, weight ~ Time | Rat)

plot(rat.fit, weight ~ Time | Diet)


###
labs <- paste(attr(BodyWeight, "labels"),
               attr(BodyWeight, "units"))
xyplot(weight ~ Time | Rat, type = c("p", "r"),
        data = BodyWeight,
        layout = c(NA, 1),  ## enforce 1 row
        groups = Diet,
        auto.key = list(title = "Diet", pch = 19, type = "n", fill = 0.6),
        pch = 19,
        xlab = labs[1],
        ylab = labs[2],
        par.settings = list(superpose.line = list(col = "black")))

# TODO
