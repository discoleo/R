########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S3:
### Composed from Simpler Subsystems
###
### draft v.0.1g


##########################
### Asymmetric Systems ###

### Composed Asymmetric Systems
### from simpler Subsystems

# - some simple Models;


###############
### History ###

### draft v.0.1f - v.0.1g:
# - system derived from Ht[3, 3, 1]:
#   [Alternating Diff] o Ht[3, 3, 1];
### draft v.0.1d - v.0.1e:
# - systems derived from Hetero-symmetric [2, 2, 1] systems:
#   x^2 + b1*y = R1;
#  -- symmetric [3, 2] sub-system; [v.0.1d]
#  -- mixt hetero-symmetric [3, 2, 1] sub-system; [v.0.1e]
### draft v.0.1b - v.0.1c:
# - first 3 sytems based on:
#   Sys 1: x^2 + y^2 + z^2 = R1; [v.0.1a - v.0.1c]
# - added more comments, structure; [v.0.1b-comm]


###################

###################
### Terminology ###
###################

### Symmetric [v, n]
# x1^n + x2^n + ... + xv^n = R1;

### Hetero-Symmetric [v, n, ...]
# x1^n + ... = R
# x2^n + ... = R
# ...
# xv^n + ... = R

#################################
#################################

#################
### Symmetric ###
#################

### Symmetric[3, 2] o Asymmetric:

### Sys 1: => {x1, x2, x3}
# x^2 + y^2 + z^2 = R1;
### Sys 2:
# x^2 + b1*x = x1
# y^2 +  x*y = x2
# z^2 +  y*z = x3

x^4 + y^4 + z^4 + 2*b[1]*x^3 + b[1]^2*x^2 + 2*x*y^3 + 2*y*z^3 + x^2*y^2 + y^2*z^2 - R[1]
x^3*y + y^3*z + b[1]*(x^2*y + x*y^2 + x*z^2) + x*y*z*(x+y+z) + (x^2*y^2 + x^2*z^2 + y^2*z^2) + b[1]*x*y*z - R[2]
b[1]*x*y*z*(x*y + x*z + y*z) + x^2*y^2*z^2 + (x^3*y^2*z + x^3*y*z^2 + x^2*y^3*z) + b[1]*x*y^3*z - R[3]

### Solution:
solve.sysEnt = function(R, b) {
	### Symmetric System: Order 2
	S = sqrt(R[1] + 2*R[2])
	S = c(S, -S)
	len = length(S)
	### x123
	x = sapply(S, function(x) roots(c(1, -x, R[2], -R[3])))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	# y, z: robust
	yz.s = S - x; yz = R[2] - x*yz.s
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2;
	z = yz.s - y;
	x = rep(x, 2); y2 = c(y, z); z = c(z, y); y = y2;
	### final subsystem
	x = sapply(as.vector(x), function(x) roots(c(1, b[1], -x)))
	y = rep(as.vector(y), each=2);
	y = sapply(seq_along(x), function(id) roots(c(1, x[id], -y[id])))
	z = rep(as.vector(z), each=4);
	z = sapply(seq_along(y), function(id) roots(c(1, y[id], -z[id])))
	x = rep(as.vector(x), each=4); y = rep(as.vector(y), each=2);
	return(cbind(x=x, y=y, z=as.vector(z)))
}

### Solution:
R = c(2, 1, -1)
b = 2

sol = solve.sysEnt(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
x^4 + y^4 + z^4 + 2*b[1]*x^3 + b[1]^2*x^2 + 2*x*y^3 + 2*y*z^3 + x^2*y^2 + y^2*z^2 # - R[1]
x^3*y + y^3*z + b[1]*x^2*y + b[1]*x*y^2 + b[1]*x*z^2 + x*y*z*(x+y+z) + x^2*y^2 + x^2*z^2 + y^2*z^2 + b[1]*x*y*z # - R[2]
b[1]*x*y*z*(x*y + x*z + y*z) + x^2*y^2*z^2 + x^3*y^2*z + x^3*y*z^2 + x^2*y^3*z + b[1]*x*y^3*z # - R[3]


#########
### Ex 2:
R = c(0, 2, -1)
b = 3

sol = solve.sysEnt(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
err = x^4 + y^4 + z^4 + 2*b[1]*x^3 + b[1]^2*x^2 + 2*x*y^3 + 2*y*z^3 + x^2*y^2 + y^2*z^2 # - R[1]
round0(err)
x^3*y + y^3*z + b[1]*x^2*y + b[1]*x*y^2 + b[1]*x*z^2 + x*y*z*(x+y+z) + x^2*y^2 + x^2*z^2 + y^2*z^2 + b[1]*x*y*z # - R[2]
x^2*y^2*z^2 + b[1]*x*y*z*(x*y + x*z + y*z) + x^3*y^2*z + x^3*y*z^2 + x^2*y^3*z + b[1]*x*y^3*z # - R[3]


### TODO:
# - evaluate techniques to simplify structure of system;


###############################

### simpler Variant:
### Sys 1: => {x1, x2, x3}
# x^2 + y^2 + z^2 = R1;
### Sys 2:
# x^2 + b1*x = x1
# y^2 +  x*y = x2
# z^2 +  x*z = x3

x^4 + y^4 + z^4 + 2*b[1]*x^3 + b[1]^2*x^2 + 2*x*y^3 + 2*x*z^3 + x^2*y^2 + x^2*z^2
x^3*y + x^3*z + x^2*y^2 + x^2*z^2 + y^2*z^2 + x*y*z*(x+y+z) + b[1]*x^2*y + b[1]*x^2*z + b[1]*x*y^2 + b[1]*x*z^2
x^2*y^2*z^2 + b[1]*x*y*z*(x*y + x*z + y*z) + x^3*y^2*z + x^3*y*z^2 + x^4*y*z + b[1]*x^3*y*z

### Pseudo-Invariants:
# s = y+z; p = y*z;
# but system can be decomposed even more extensively;
x^4 + y^4 + z^4 + 2*b[1]*x^3 + b[1]^2*x^2 + 2*x*(s^3 - 3*s*p) + x^2*(s^2 - 2*p)
x^3*s + x^2*(s^2 - p + b[1]*s) + p^2 + x*p*s + b[1]*x*(s^2 - 2*p)
x^4*p + x^3*p*s + x^2*p^2 + b[1]*x^3*p + b[1]*x^2*p*s + b[1]*x*p^2
(x^2 + s*x + p)*(x^2*p + b[1]*x*p) # rewritten Eq 3

### Solution:
solve.sysEnt = function(R, b) {
	### Symmetric System: Order 2
	S = sqrt(R[1] + 2*R[2])
	S = c(S, -S)
	len = length(S)
	### x123
	x = sapply(S, function(x) roots(c(1, -x, R[2], -R[3])))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	# y, z: robust
	yz.s = S - x; yz = R[2] - x*yz.s
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2;
	z = yz.s - y;
	x = rep(x, 2); y2 = c(y, z); z = c(z, y); y = y2;
	### final subsystem
	x = sapply(as.vector(x), function(x) roots(c(1, b[1], -x)))
	y = rep(as.vector(y), each=2);
	y = sapply(seq_along(x), function(id) roots(c(1, x[id], -y[id])))
	z = rep(as.vector(z), each=2);
	z = sapply(seq_along(x), function(id) roots(c(1, x[id], -z[id])))
	# TODO: all (y, z) combinations
	x = rep(as.vector(x), each=2); # y = rep(as.vector(y), each=2);
	return(cbind(x=x, y=as.vector(y), z=as.vector(z)))
}
test = function(x, y, z, b, R) {
	err1 = x^4 + y^4 + z^4 + 2*b[1]*x^3 + b[1]^2*x^2 + 2*x*y^3 + 2*x*z^3 + x^2*y^2 + x^2*z^2
	err2 = x^3*y + x^3*z + x^2*y^2 + x^2*z^2 + y^2*z^2 + x*y*z*(x+y+z) + b[1]*x^2*y + b[1]*x^2*z + b[1]*x*y^2 + b[1]*x*z^2
	err3 = x^2*y^2*z^2 + b[1]*x*y*z*(x*y + x*z + y*z) + x^3*y^2*z + x^3*y*z^2 + x^4*y*z + b[1]*x^3*y*z
	err = rbind(err1, err2, err3)
	if( ! missing(R)) {
		# err = sapply(seq(nrow(err)), function(id) err[id,] - R[id])
		# err = t(err) # which is better?
		err = sapply(seq(ncol(err)), function(id) err[,id] - R)
	}
	round0(err)
}

### Solution:
R = c(2, 1, -1)
b = 2

sol = solve.sysEnt(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
test(x, y, z, b=b)

poly.calc(x[seq(1,47, by=2)])


###############################
###############################

### Symmetric[3, 2] o Assymetric:
### Sys 1: => {x1, x2, x3}
# x^2 + y^2 + z^2 = R1;
### Sys 2:
# x+y+z = x1
# x-y+z = x2
# x+y-z = x3

### Entangled system:
3*x^2 + 3*y^2 + 3*z^2 + 2*x*y + 2*x*z - 2*y*z - R[1]
3*x^2 - y^2 - z^2 + 2*x*y + 2*x*z + 2*y*z - R[2]
x^3 - y^3 - z^3 + x^2*y + x^2*z - x*y^2 + y^2*z + y*z^2 - x*z^2 + 2*x*y*z - R[3]

### Diff: [1] - [2]
y^2 + z^2 - y*z - (R[1]-R[2])/4 # = 0
# Eq 2 =>
3*x^2 + 2*x*y + 2*x*z + y*z - (R[1] + 3*R[2])/4 # = 0
# Eq 3 =>
x^3 + x^2*y + x^2*z - x*y^2 + y^2*z + y*z^2 - x*z^2 + 2*x*y*z - (y+z)*(R[1] - R[2])/4 - R[3] # = 0

### TODO:
# - transform system into more compact equations;

### Solution:
solve.sysEnt = function(R, b) {
	### Symmetric System: Order 2
	S = sqrt(R[1] + 2*R[2])
	S = c(S, -S)
	len = length(S)
	### x123
	x = sapply(S, function(x) roots(c(1, -x, R[2], -R[3])))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	# y, z: robust
	yz.s = S - x; yz = R[2] - x*yz.s
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2;
	z = yz.s - y;
	x = as.vector(rep(x, 2)); y2 = c(y, z); z = as.vector(c(z, y)); y = as.vector(y2);
	### final subsystem
	
	return(cbind(x=(y+z)/2, y=(x-y)/2, z=(x-z)/2))
}
test = function(x, y, z, R) {
	err1 = 3*x^2 + 3*y^2 + 3*z^2 + 2*x*y + 2*x*z - 2*y*z
	err2 = 3*x^2 - y^2 - z^2 + 2*x*y + 2*x*z + 2*y*z
	err3 = x^3 - y^3 - z^3 + x^2*y + x^2*z - x*y^2 + y^2*z + y*z^2 - x*z^2 + 2*x*y*z
	err = rbind(err1, err2, err3)
	if( ! missing(R)) {
		# err = sapply(seq(nrow(err)), function(id) err[id,] - R[id])
		# err = t(err) # which is better?
		err = sapply(seq(ncol(err)), function(id) err[,id] - R)
	}
	round0(err)
}

### Solution:
R = c(2, 1, -1)

sol = solve.sysEnt(R)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
test(x, y, z)



#########
### TODO:
x^4*z^2 + x^2*y^4 + y^2*z^4 + b1*x*y + b1*x*z + b1*y*z + 2*x*y^3*z^2 + 2*x^2*y*z^3 + 2*x^3*y^2*z
x^2*y^2 + x^2*z^2 + y^2*z^2 + 2*x^2*y*z + 2*x*y^2*z + 2*x*y*z^2 + b1*x*y*z
b1*y*z^2 + b1*x*y^2 + b1*x^2*z + x^2*y^2*z^2


###############################
###############################

##################
### Base: Ht 2 ###
##################

### Hetero-Symmetric[2, 2, 1] o
### o Symetric[3, 2]:
### Sys 1: => {xi, yi}
# x^2 + b1*y = R1;
# y^2 + b1*x = R1;
### Sys 2:
# x^2 + y^2 + z^2 = xi
# x*y + x*z + y*z = yi
# x*y*z = R2

### symmetric System:
x^4 + y^4 + z^4 + b[1]*(x*y + x*z + y*z) + 2*x^2*y^2 + 2*x^2*z^2 + 2*y^2*z^2 - R[1]
b[1]*(x^2 + y^2 + z^2) + 2*x*y*z*(x+y+z) + x^2*y^2 + x^2*z^2 + y^2*z^2 - R[1]
x*y*z - R[2]
# system is extensively decomposable into Sys1 o Sys2;

### Classic:
### Eq 1:
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 + b[1]*E2 + 2*(E2^2 - 2*E3*S) - R[1]
S^4 - 4*E2*S^2 + 4*E2^2 + b[1]*E2 - R[1]
### Eq 2:
b[1]*(S^2 - 2*E2) + 2*E3*S + (E2^2 - 2*E3*S) - R[1]
b[1]*S^2 + E2^2 - 2*b[1]*E2 - R[1]
# =>
b[1]^2*S^4 - 4*b[1]^2*E2*S^2 + 4*b[1]^2*E2^2 + b[1]^3*E2 - b[1]^2*R[1]
(E2^2 - 2*b[1]*E2 - R[1])^2 + 4*b[1]*E2*(E2^2 - 2*b[1]*E2 - R[1]) + 4*b[1]^2*E2^2 + b[1]^3*E2 - b[1]^2*R[1]
E2^4 - 2*R[1]*E2^2 + b[1]^3*E2 + R[1]^2 - b[1]^2*R[1]
(E2^2 - R[1])^2 + b[1]^2*(b[1]*E2 - R[1])
# - both the equations for E2 & S are order 4 polynomials;
# - the system can be actually decomposed further;


### Solution:
solve.sysEntHt2 = function(R, b) {
	### Hetero-Symmetric System: Order 2
	sol1 = roots(c(1, b[1], -R[1])) # roots x == y;
	S = b[1] # 2nd set of roots (x != y);
	xy = - R[1] + (S^2 + b[1]*S)/2
	xy.d = sqrt(S^2 - 4*xy + 0i)
	x = (S + xy.d) / 2; x = c(x, S-x); y = S - x;
	RS = cbind(x, y); RS = rbind(RS, cbind(sol1, sol1))
	### Symmetric System: Order 2
	S = sqrt(RS[,1] + 2*RS[,2])
	S = c(S, -S); RS = rbind(RS, RS)
	len = length(S)
	### x123
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], RS[id,2], -R[2])))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	RS = apply(RS, 2, function(x) rep(x, each=3))
	# y, z: robust
	yz.s = S - x; yz = RS[,2] - x*yz.s
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2;
	z = yz.s - y;
	x = as.vector(rep(x, 2)); y2 = c(y, z); z = as.vector(c(z, y)); y = as.vector(y2);
	
	return(cbind(x=x, y=y, z=z))
}
test = function(x, y, z, b, R) {
	err1 = x^4 + y^4 + z^4 + b[1]*(x*y + x*z + y*z) + 2*x^2*y^2 + 2*x^2*z^2 + 2*y^2*z^2
	err2 = b[1]*(x^2 + y^2 + z^2) + 2*x*y*z*(x+y+z) + x^2*y^2 + x^2*z^2 + y^2*z^2
	err3 = x*y*z
	err = rbind(err1, err2, err3)
	if( ! missing(R)) {
		R = c(R[1], R[1], R[length(R)])
		# err = sapply(seq(nrow(err)), function(id) err[id,] - R[id])
		# err = t(err) # which is better?
		err = sapply(seq(ncol(err)), function(id) err[,id] - R)
	}
	round0(err)
}

### Solution:
R = c(-2, 1)
b = 3

sol = solve.sysEntHt2(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
test(x, y, z, b)

# polynomial of order 48
round0.p(poly.calc(x))
# same roots for initial system (25:48) can be factored
round0.p(poly.calc(x[1:24]))


###############################
###############################

### Hetero-Symmetric[2, 2, 1] o
### o Mixt-HeteroSymetric[3, 2, 1]

### Sys 1: => {xi, yi}
# x^2 + b1*y = R1;
# y^2 + b1*x = R1;
### Sys 2:
# x*y^2 + y*z^2 + z*x^2 = xi
# x*y + x*z + y*z = yi
# x*y*z = R2

x^4*z^2 + x^2*y^4 + y^2*z^4 + b[1]*(x*y + x*z + y*z) + 2*x^3*y^2*z + 2*x*y^3*z^2 + 2*x^2*y*z^3 - R[1]
b[1]*(x^2*z + x*y^2 + y*z^2) + 2*x*y*z*(x + y + z) + (x^2*y^2 + x^2*z^2 + y^2*z^2) - R[1]
x*y*z - R[2]


### Solution:
solve.ht3 = function(R, b=0) {
	if(length(b) == 1 && b[1] == 0) {
		coeff = c(R[3], 0, - (R[1]+6*R[3])*R[2], R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
	} else {
		coeff = c(R[3], (b[1]^2 + b[1]*R[2]), - (R[1]*R[2] + 3*b[1]*R[3] + 6*R[2]*R[3] + 2*b[1]*R[1]),
			R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
	}
	if(length(b) > 1) {
		# Ext 2:
		coeff = coeff - c(b[2]^3 + b[1]*b[2], -b[2]*(R[1] + 3*b[2]*R[2] + 6*R[3]), 3*b[2]*R[2]^2, 0)
	}
	S = roots(coeff)
	b2 = if(length(b) > 1) b[2] else 0; # Ext 2;
	x = sapply(S, function(x) roots(c(1,-x, R[2] - b2*x, -R[3])))
	S = matrix(S, ncol=3, nrow=3, byrow=T)
	yz = R[3]/x
	yz.s = S - x
	### robust:
	y = - (x*yz - (x^2+yz)*yz.s - b[1]*S + R[1]) / (x^2+yz - x*yz.s)
	z = yz.s - y
	sol = cbind(as.vector(x), as.vector(y), as.vector(z))
	return(sol)
}
solve.sysEnt = function(R, b) {
	### Hetero-Symmetric System: Order 2
	sol1 = roots(c(1, b[1], -R[1])) # Set 1: roots x == y;
	S = b[1] # Set 2: roots (x != y);
	xy = - R[1] + (S^2 + b[1]*S)/2
	xy.d = sqrt(S^2 - 4*xy + 0i)
	x = (S + xy.d) / 2; x = c(x, S-x); y = S - x;
	RS = cbind(x, y); RS = rbind(RS, cbind(sol1, sol1))
	### Mixt Hetero-Symmetric System: Order 2
	sol = lapply(seq(nrow(RS)), function(id) solve.ht3(c(RS[id,1], RS[id,2], R[2]), b=0))
	sol = do.call(rbind, sol)
	return(sol)
}
test = function(x, y, z, b, R) {
	err1 = x^4*z^2 + x^2*y^4 + y^2*z^4 + b[1]*(x*y + x*z + y*z) + 2*x^3*y^2*z + 2*x*y^3*z^2 + 2*x^2*y*z^3
	err2 = b[1]*(x^2*z + x*y^2 + y*z^2) + 2*x*y*z*(x + y + z) + (x^2*y^2 + x^2*z^2 + y^2*z^2)
	err3 = x*y*z
	err = rbind(err1, err2, err3)
	if( ! missing(R)) {
		R = c(R[1], R[1], R[length(R)])
		# err = sapply(seq(nrow(err)), function(id) err[id,] - R[id])
		# err = t(err) # which is better?
		err = sapply(seq(ncol(err)), function(id) err[,id] - R)
	}
	round0(err)
}

### Example:
R = c(-2, 1)
b = 1

sol = solve.sysEnt(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
test(x, y, z, b)

round0.p(poly.calc(x[1:18]))
# the 2nd set of 2*9 roots can be factored;

err = 1 - 3*x + 12*x^2 - 25*x^3 + 39*x^4 + 15*x^5 - 26*x^6 + 96*x^7 + 114*x^8 + 80*x^9 + 186*x^10 +
	+ 198*x^11 + 143*x^12 + 96*x^13 + 57*x^14 + 2*x^15 - 9*x^16 + x^18
round0(err) # only 18 are still roots
# but there is also some numerical inaccuracy!


### Ex 2:
R = c(0, 1)
b = 3

sol = solve.sysEnt(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
test(x, y, z, b)

round0.p(poly.calc(x[1:18]))
# the 2nd set of 2*9 roots can be factored;

err = 1 - 9*x + 54*x^2 - 195*x^3 + 495*x^4 - 567*x^5 + 132*x^6 + 1350*x^7 + 1006*x^9 + 2466*x^10 +
	+ 1242*x^11 + 1365*x^12 + 828*x^13 + 243*x^14 - 42*x^15 - 27*x^16 + x^18
round0(err) # only 18 are still roots


###############################
###############################

##################
### Base: Ht 3 ###
##################

### Order 3

### Base System: Ht[3,3,1]
### System 2: Ht-Diff
# Ht-Diff o Ht[3,3,1]
# (x,y,z) => - x + y + z, x - y + z, x + y - z;
# Transform: Eq[1] + Eq[2], ...;

x^3 - 6*x*y*z + 3*x*y^2 + 3*x*z^2 + b1*y = R
y^3 - 6*x*y*z + 3*y*z^2 + 3*x^2*y + b1*z = R
z^3 - 6*x*y*z + 3*y^2*z + 3*x^2*z + b1*x = R


### TODO: (-2,2,1)
193*x^3 - 4632*x*y*z + 42*b1*x + 97*b1*y + 54*b1*z + 1818*x^2*y - 666*x^2*z - 132*x*y^2 + 2172*x*z^2 + 1800*y^2*z - 360*y*z^2
193*y^3 - 4632*x*y*z + 54*b1*x + 42*b1*y + 97*b1*z + 2172*x^2*y - 360*x^2*z - 666*x*y^2 + 1800*x*z^2 + 1818*y^2*z - 132*y*z^2
193*z^3 - 4632*x*y*z + 97*b1*x + 54*b1*y + 42*b1*z + 1800*x^2*y - 132*x^2*z - 360*x*y^2 + 1818*x*z^2 + 2172*y^2*z - 666*y*z^2

### TODO: (-3,1,1)
700*x^3 + 504*x*y*z - 24*b1*x + 76*b1*y - 24*b1*z - 696*x^2*y - 696*x^2*z + 204*x*y^2 + 204*x*z^2 - 96*y^2*z - 96*y*z^2
700*y^3 + 504*x*y*z - 24*b1*x - 24*b1*y + 76*b1*z + 204*x^2*y - 96*x^2*z - 696*x*y^2 - 96*x*z^2 - 696*y^2*z + 204*y*z^2
700*z^3 + 504*x*y*z + 76*b1*x - 24*b1*y - 24*b1*z - 96*x^2*y + 204*x^2*z - 96*x*y^2 - 696*x*z^2 + 204*y^2*z - 696*y*z^2




###############################
###############################

### TODO:
# x+y, y+z, z-x
x^2 + y^2 + z^2 + x*y - x*z + y*z
- x^2 + y^2 + z^2 - x*y + x*z + 3*y*z
x*z^2 - x^2*y - x^2*z - x*y^2 + y^2*z + y*z^2


### TODO:
# HtAs[3, 2, 1] & x+y, y+z, z-x: Paradox
x^2 + y^2 + b[1]*y + b[1]*z + 2*x*y
y^2 + z^2 - b[1]*x + b[1]*z + 2*y*z
x^2 + z^2 + b[1]*x + b[1]*y - 2*x*z
# =>
x^2 + y^2 + z^2 + b1*y + b1*z + x*y - x*z + y*z = 3/2*R

### Solution:
# TODO: Paradox
solve.sysHt32 = function(R, b) {
	x.sum = roots(c(2, b[1], (b[1]^2 - 2*R[1]), - 3*b[1]*R[1] + 6*b[1]^3))
	E2 = (x.sum^2 + b[1]*x.sum - 3*R[1])/2
	E3 = E2*x.sum + b[1]^3
	x = sapply(1:length(x.sum), function(id) roots(c(1, -x.sum[id], E2[id], -E3[id])))
	x = cbind(as.vector(x[,-1])) # TODO: remove robustly the set of wrong solutions
	y = (R[1] - x^2)/b[1]
	z = (R[1] - y^2)/b[1]
	sol = cbind(x=x, y=y, z=z)
	sol
}
solve.sysEntHtLSD = function(R, b) {
	# uses solve.sysHt32()!
	sol = solve.sysHt32(R, b)
	x = 1
	y = sol[,1] - x
	z = sol[,3] + x
	cbind(x=x, y=y, z=z, sol)
}

### Example:
R = 1
b = 3

sol = solve.sysEntHtLSD(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

###############################
###############################

####################
### Experimental ###

### TODO:
# - explore, clean;

x = sqrt(2); y = sqrt(3); z = sqrt(5);
E2 = x*y + x*z + y*z;
E3 = x*y*z;
S = sum(x,y,z)
#
b = 1
n = 2
# R[i]
R1 = x^n - b*y
R2 = y^n - b*z
R3 = z^n - b*x
R12 = R1*R2 + R1*R3 + R2*R3

### Test: n = 2
S^2 - 2*E2 - b*S - (R1+R2+R3) # = 0

(x^3+y^3+z^3) - b*E2 - (R1*x + R2*y + R3*z)
S^3 - 3*E2*S + 3*E3 - b*E2 - (R1*x + R2*y + R3*z) # = 0

(E2^2 - 2*E3*S) - b*(S^3 - 3*E2*S + 3*E3) - b*(R1*z + R2*x + R3*y) - R12

(x^2*y^2 + x^2*z^2 + y^2*z^2) - b*(x^3+y^3+z^3) - (R1*y^2 + R2*z^2 + R3*x^2)
(E2^2 - 2*E3*S) - b*(S^3 - 3*E2*S + 3*E3) - (R1*y^2 + R2*z^2 + R3*x^2)


#############

### x+y
2*(x*y^2 + y*z^2 + x^2*z) + 3*(y^2*z + x*z^2 + x^2*y) + x^3 + y^3 + z^3 + 6*x*y*z
3*y*z + 3*x*y + 3*x*z + x^2 + y^2 + z^2
y*z^2 + y^2*z + x*y^2 + x*z^2 + x^2*y + x^2*z + 2*x*y*z

### =>
x^3 + y^3 + z^3 - (x*y^2 + y*z^2 + x^2*z) - (R1 - 3*R3)
x^2 + y^2 + z^2 + 3*E2 - R2


### TODO:
(x*y^2 + y*z^2 + x^2*z) + b*(y^2*z + x*z^2 + x^2*y) = R1



