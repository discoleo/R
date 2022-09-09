
### Exam Questions for putative students

### Activation Functions

# Q1.) Explain the characteristics of the following (activation) functions:

lim = c(-3, 3)
curve((tanh(x) + 1)/2, lim[1], lim[2], ylim=c(-0.125,1.125))
curve((tanh(x) + tanh(x+1))^2/4, add=TRUE, col="red")
curve((tanh(x-1) + tanh(x+1))^2/4, add=TRUE, col="blue")
curve((tanh(x-2) + tanh(x+1))^2/4, add=TRUE, col="#D09632")

abline(h=1, col="grey")


### Probabilities

# Q2.) Explain the concept of probability density function (PDF) and the normalization of the PDF.
#  Can any of these functions represent a PDF?

# Q3.) Explain the concept of cumulative distribution function (CDF).
#  Can any of these functions represent a CDF?

# Q4.) Create the functions G[i] = 1 - f[i](x), where f[i] represents
#  the last 3 functions defined above ( i = 1 to 3). Compare the G[i] functions
#  to a set of Gaussian distributions. Explain if these functions behave similarly
#  to a Gaussian distribution or if there are relevant differences?

