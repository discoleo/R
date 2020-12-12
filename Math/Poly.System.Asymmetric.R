########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S3:
### Composed from Simpler Subsystems
###
### draft v.0.1c


##########################
### Asymmetric Systems ###

### Composed Asymmetric Systems
### from simpler Subsystems


###############
### History ###

### draft v.0.1b - v.0.1c:
# - first 3 sytems based on:
#   Sys 1: x^2 + y^2 + z^2 = R1; [v.0.1a - v.0.1c]
# - added more comments, structure; [v.0.1b-comm]


#################################

### Symmetric[3, 2] o Assymetric:
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



