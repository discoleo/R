#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Log( F(Y) )-Types


###################
### Log-Product ###
###################

### log(y + x) * log(y - x) = P(x)

### D =>
(dy + 1)*(y - x)*log(y - x) + (dy - 1)*(y + x)*log(y + x) - (y^2 - x^2)*dp # = 0

### D2 =>
(d2y*(y - x) + dy^2 - 1)*log(y - x) + (d2y*(y + x) + dy^2 - 1)*log(y + x) +
	+ 2*dy^2 - 2 - 2*(y*dy - 2*x)*dp - (y^2 - x^2)*d2p # = 0

### Solve linear system:
log(y + x) = - ... /2 / (y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)
log(y - x) = ... /2 / (y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)


### ODE:
# TODO

### Example:
### P(x) = x
(dy + 1)*(y - x)*log(y - x) + (dy - 1)*(y + x)*log(y + x) - (y^2 - x^2) # = 0
(d2y*(y - x) + dy^2 - 1)*log(y - x) + (d2y*(y + x) + dy^2 - 1)*log(y + x) +
	+ 2*(dy^2 - y*dy + x - 1) # = 0
# =>
2*log(y + x) = - (y - x) * (2*(dy + 1)*(dy^2 - y*dy + x - 1) + (d2y*(y - x) + dy^2 - 1)*(y + x)) /
	(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)
2*log(y - x) = - (y + x) * (2*(dy - 1)*(dy^2 - y*dy + x - 1) + (d2y*(y + x) + dy^2 - 1)*(y - x)) /
	(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)
# =>
2*log(y + x) = - (y - x)*(y^2*d2y - x^2*d2y + 2*dy^3 - y*dy^2 + x*dy^2 + 2*dy^2 - 2*y*dy + 2*x*dy - 2*dy - y + x - 2) /
	(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)
2*log(y - x) = - (y + x)*(...) /
	(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)

### ODE:
(y^2 - x^2)*(...)*(...) +
	- 4*x*(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)^2 # = 0

### TODO: check!

