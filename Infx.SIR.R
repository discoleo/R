
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
plot.sir = function(y, times, legend.xy=c(max(times)*2/3, 0.7)) {
	matplot(x = times, y = y, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
        lwd = 1, lty = 1, bty = "l", col = 2:6)

	## Add legend
	legend(legend.xy[1], legend.xy[2], c("Susceptible", "Infected", "Recovered", "Terminal", "Dead"),
		pch = 1, col = 2:6, bty = "n")
}
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

    dS = -beta * S * I - sv * S * T - vi * S * T^2 - vc * S * T * D
    dI =  beta * S * I - gamma * I - m * I + sv * S * T + vi * S * T^2 + vc * S * T * D 
    dR =  gamma * I
	dT = m * I - mt * T
	dD = mt * T

    return(list(c(dS, dI, dR, dT, dD)))
  })
}

### Set parameters

## Proportion in each compartment:
# Susceptible = 0.999999;
# Infected = 0.000001;
# Recovered = 0;
init = c(S = 1-1e-6, I = 1e-6, R = 0.0, T = 0.0, D = 0.0)

## beta: infection parameter;
## gamma: recovery parameter;
beta = 1.4247 / 4
gamma = 0.14286
#
parameters = c(beta = beta, gamma = gamma, sv = 3*beta,
	vi = 2, vc = 10,
	m = gamma * 10000/80000, mt = gamma)
## Time frame
end.time = 120 # 70
times = seq(0, end.time, by = 1)


## Solve using ode
out = solve.sir(sir, init, parameters, times)
head(out, 10)


#############

# png(file="SIR.AdvancedModel.png", bg ="white", pointsize=15)

plot.sir(out, times, legend.xy=c(end.time*2/3, 0.7))
# essential function for any scientific publishing
plot.leo(end.time=end.time)

# dev.off()

