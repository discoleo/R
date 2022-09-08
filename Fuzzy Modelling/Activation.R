

### Activation Functions

# 1.) Explain the characteristics of the following (activation) functions:

lim = c(-3, 3)
curve((tanh(x) + 1)/2, lim[1], lim[2], ylim=c(-0.125,1.125))
curve((tanh(x) + tanh(x+1))^2/4, add=TRUE, col="red")
curve((tanh(x-1) + tanh(x+1))^2/4, add=TRUE, col="blue")
curve((tanh(x-2) + tanh(x+1))^2/4, add=TRUE, col="#D09632")

abline(h=1, col="grey")

