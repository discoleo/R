

### Linearisation vs NLS

### Example: Michaelis-Menten Equation;
# V = Vmax * s / (Km + s)


### Discussion:

### Issues with Linearisation:
# The Michaelis-Menten equation can be easily linearised.
# However, the linearisation creates various other problems:
# - Non-Linearized eq differs substantially for larger values of s;
# - Transformed values 1/s and 1/V are relatively small when "s" is large;
# - Linearization fails therefore for these upper values;

### Workaround:
# - it is possible to use weights during the fitting process;
# - Weights = s virtually eliminates the effect of linearisation;
# - However, the larger values of s may be noisy as well;
# - Weights = sqrt(s) is probably more robust;

####################

### Helper Functions

library(ggplot2)

### Linearisation:
# w = weights;
predict.enz = function(eq, data = enz.data, print = TRUE, w = NULL) {
	if(is.null(w)) {
		enz.fit = lm(eq, data=data);
	} else {
		data$w = w;
		enz.fit = lm(eq, data=data, weights = w);
	}
	if(print) print(summary(enz.fit));
	fit.data = predict(enz.fit, data);
	model.dt = cbind(data, .fitted = fit.data);
	# model.dt = broom::augment(enz.fit, data=data);
	# names(model.dt)[1] = "s";
	invisible(model.dt);
}

### Plot:
# FUN = function to transform fitted data;
# ... = additional parameters to geom_line;
plot.enz = function(data, gg = NULL, FUN = function(x) x,
		fitted = TRUE, col.fit = "black", size = 1.4, size.p = 2,
		lbl  = c("s (mM)", "V (mM / s)"), ...) {
	names(data)[1:2] = c("x", "y");
	if(is.null(gg)) {
		gg = ggplot(data = data) +
			geom_point(aes(x = x, y = y),
				color = "red", size = size.p) +
			labs(
				x = lbl[1], y = lbl[2]);
	}
	if(fitted & ".fitted" %in% names(data)) {
		gg = gg + geom_line(
			aes(x = x, y = FUN(.fitted)),
			data = data, color = col.fit, linewidth = size, ...);
	}
	return(gg);
}


####################

### Ref:
# Based on: J. Zobitz [2023]
# Exploring Modeling with Data and Differential Equations Using R
# ISBN 9781032259482
# Ch 8: Linear Regression and Curve Fitting
# https://jmzobitz.github.io/ModelingWithR/linear-regression-08.html

### Enzyme Kinetic:
enz.data = data.frame(
	s = c(0.1, 0.2, 0.5, 1.0, 2.0, 3.5, 5.0),
	V = c(0.04, 0.08, 0.17, 0.24, 0.32, 0.39, 0.42)
)

### Linear Model:
enz.eq = I(1 / V) ~ 1 + I(1 / s)

### NLS Model:
nls.eq = V ~ Vm * s / (Km + s);


####################
# Testing Linearity:
# - Reciprocals: reasonably linear;
ggplot(data = enz.data) +
	geom_point(aes(x = 1 / s, y = 1 / V),
		color = "red", size = 2) +
	labs(
		x = "1/s (1/mM)",
		y = "1/V (s / mM)")


### Linear Model:
enz.fit = lm(enz.eq, data = enz.data)
summary(enz.fit)


# Add fitted data to the model
# - using broom::augment:
# enz.dm = broom::augment(enz.fit, data = enz.data);
# - using native R:
enz.dm = cbind(enz.data, .fitted = predict(enz.fit, enz.data));

### Plot:
# Compare fitted model to the *linearised* data:
ggplot(data = enz.data) +
	geom_point(aes(x = 1 / s, y = 1 / V),
		color = "red", size = 2) +
	geom_line(
		data = enz.dm,
		aes(x = 1 / s, y = .fitted)) +
	labs(
		x = "1/s (1/mM)",
		y = "1/V (s / mM)")


### Native Data:
# [Shorthand]
enz.dm = predict.enz(enz.eq, data = enz.data);
plot.enz(enz.dm, FUN = \(x) 1/x)

# Note:
# - Linearization fails for the larger values of s;


#####################

### Weighted Approach

col = c("#FA44B2A0", "#E29632A0", "#32FF44A0")

# Linear:
enz.dm  = predict.enz(enz.eq, data = enz.data);
enz.dmw = predict.enz(enz.eq, data = enz.data, w = sqrt(enz.data$s))
enz.dmx = predict.enz(enz.eq, data = enz.data, w = enz.data$s)
# NLS:
enz.nls = nls(nls.eq,
	data  = enz.data, control = nls.control(tol = 1E-8, minFactor = 1/2^14),
	start = list(Km = 1.2, Vm = 0.533));
enz.dnl = cbind(enz.data, .fitted = predict(enz.nls, enz.data))

# Plot:
gg = plot.enz(enz.dm, FUN = \(x) 1/x); # Non-Wieghted
gg = plot.enz(enz.dmw, gg=gg, FUN = \(x) 1/x, col = col[1])
gg = plot.enz(enz.dmx, gg=gg, FUN = \(x) 1/x, col = col[2])
gg = plot.enz(enz.dnl, gg=gg, col = col[3])
gg;


# Note:
# - Weights = s fits almost perfectly;
# - BUT: is probably far more sensitive to noisy data
#   (e.g. when large "s" is noisy);

