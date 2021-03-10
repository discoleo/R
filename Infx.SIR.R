
#############################
###
### SIR Model
### Extension of the model
###
### Leonard Mada

# Features:
# - added mortality:
#   T = terminal patients, D = dead;
# - m = mortality, mt = time to death;
#   [but NO birth-rate]
# - added higher-order interaction terms:
#  -- sv = super-virulence of terminal patients;
#  -- vi = interactions between virulence factors;
#  -- vc = cummulative virulence,
#     e.g. due to cummulative mutations acquired by the disease;

# - code based on:
#   https://archives.aidanfindlater.com/blog/2010/04/20/the-basic-sir-model-in-r/
#   http://rstudio-pubs-static.s3.amazonaws.com/6852_c59c5a2e8ea3456abbeb017185de603e.html
# - model extended by LM;

##########

# install.packages("deSolve")

library(deSolve)


### Solve SIR
solve.sir = function(sir.f, init, parameters, times) {
	## Solve using ode (General Solver for Ordinary Differential Equations)
	out <- ode(y = init, times = times, func = sir.f, parms = parameters)
	## change to data frame
	out <- as.data.frame(out)
	## Delete time variable
	out$time <- NULL
	return(out)
}

### Plot SIR
basic.lbl = c("Susceptible", "Infected", "Recovered");
legend.xyf = function(times, y=0.7, off=c(0,0)) {
	c(max(times)*2/3, y) + off;
}
plot.sir = function(y, times, legend.lbl=c(basic.lbl, "Terminal", "Dead"),
		legend.xy, leg.off=c(0,0), ylab="Susceptible and Recovered", col) {
	if(missing(legend.xy)) legend.xy=legend.xyf(times, y=max(y) * 0.7, leg.off)
	if(missing(col)) col = seq(2, ncol(y) + 1);
	matplot(x = times, y = y, type = "l",
        xlab = "Time", ylab = ylab, main = "SIR Model",
        lwd = 1, lty = 1, bty = "l", col = col)

	## Add legend
	legend(legend.xy[1], legend.xy[2], legend.lbl,
		pch = 1, col = col, bty = "n")
}
# critical ploting function
plot.leo = function(x=NA, y=0.25, end.time=NA, adj=0, addYear=TRUE) {
	# critical function for anyone wanting to publish in Nature
	str = "Created\nby Leo"
	if(addYear) {
		str = paste(str, " (", format(Sys.Date(), "%Y"), ")", sep="", collapse="")
	}
	if( is.na(x) && ! is.na(end.time) ) {
		x = end.time*2/3 + 9
	}
	text(x, y, str, adj=adj)
}

### Create an SIR function
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
	
	dVirulence = vi * S * T^2; # has almost NO effect in this model;
	dcVir = vc * S * I * D;
    dS = -beta * S * I - sv * S * T - dVirulence - dcVir;
    dI =  beta * S * I - gamma * I - m * I + sv * S * T + dVirulence + dcVir; 
    dR =  gamma * I
	dT = m * I - mt * T
	dD = mt * T

    return(list(c(dS, dI, dR, dT, dD)))
  })
}

### Set parameters

### Proportion in each compartment:
# Susceptible = 0.999999;
# Infected = 0.000001;
# Recovered = 0;
init = c(S = 1-1e-6, I = 1e-6, R = 0.0, T = 0.0, D = 0.0)


### Time frame
end.time = 120 # 70
times = seq(0, end.time, by = 1)

### Parameters
# beta: infection parameter;
# gamma: recovery parameter;
beta = 1.4247 / 4
gamma = 0.14286 / 1.2
#
parameters = c(beta = beta, gamma = gamma, sv = 3*beta,
	vi = 2, vc = 10,
	m = gamma * 10000/80000, mt = gamma/2)


## Solve using ode
out = solve.sir(sir, init, parameters, times)
head(out, 10)


### Plot

# png(file="SIR.AdvancedModel.png", bg ="white", pointsize=15)

plot.sir(out, times, legend.xy=c(end.time*2/3, 0.7))
# essential function for any scientific publishing
plot.leo(end.time=end.time)

# dev.off()

###################
###################

### Comparison

solve.SIR = function(param, parameters, type="sv") {
	type = match(type, c("beta", "gamma", "sv", "vi", "vc", "m", "mt"))
	parameters[type] = param;
	solve.sir(sir, init, parameters, times)$D;
}

sv.seq = seq(0, 1, by=0.2)
out = sapply(sv.seq, solve.SIR, parameters=parameters, type="sv")
plot.sir(out, times, ylab="Death", legend.lbl=paste("SV", sv.seq))


# very small effect: dying earlier;
vc.seq = seq(0, 5000, by=1000)
out = sapply(vc.seq, solve.SIR, parameters=parameters, type="vc")
plot.sir(out, times, ylab="Death", legend.lbl=paste("CummVir", vc.seq))

