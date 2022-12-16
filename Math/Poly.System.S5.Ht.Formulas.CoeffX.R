########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric / Formulas
### X-Coeffs
###
### draft v.0.1a


### X-Coefficients:
# - used for computing x2 & x3;
# - System: Mixed S5Ht, Base system;


# this file:
# source("Poly.System.S5.Ht.Formulas.CoeffX.R")


# e = matrix with the e-components;
# - used in: Poly.System.S5.Ht.Formulas.R;
coeff.S5HtMixed.x2 = function(x1, e, E11a, E11b) {
	s = e[,1]; e2 = e[,2]; e3 = e[,3]; e4 = e[,4];
	### Coeff b:
	b20r = s^2 - 2*s*x1 + x1^2 + 2*E11b - e2;
	b02r = x1^2 + E11a + 2*E11b - e2;
	b10r = s^2*x1 - s*x1^2 - 2*s*E11b + x1*E11b + e3;
	b01r = s^2*x1 + 2*s*x1^2 - 2*E11a*x1 - s*E11b - 3*x1*E11b + e3;
	b00r = - s*x1^3 + E11a*x1^2 + x1^2*E11b - s*x1*E11b + E11b^2 - 2*e4;
	b11r = s^2 + 3*E11b;
	### Coeff c:
	c00 = e3*x1^2;
	c10 = E11b^2 - s*E11b*x1;
	c01 = - (s*E11b*x1 - E11b^2 + 2*e3*x1);
	c20 = s^2*x1 - 2*s*E11b + E11b*x1;
	c02 = s^2*x1 + s*x1^2 - s*E11b - E11b*x1 + e3;
	c11 = 2*s^2*x1 - 3*s*E11b;
	c21 = 2*s^2 - 4*s*x1 + 5*E11b;
	c12 = s^2 - s*x1 + 4*E11b;
	c30 = s^2 - 2*s*x1 + 2*E11b;
	c03 = E11b - s*x1;
	#
	c00r = c00 - e4*(x1 - s);
	c10r = c10 - 3*e4 - e3*s + e3*x1;
	c01r = c01 - 4*e4;
	c11r = c11 + 6*e3; c02r = c02;
	c20r = c20 + e3 + e2*s - e2*x1;
	c21r = c21 - 4*e2; c12r = c12 - 2*e2;
	c30r = c30 - s^2 - e2 + x1*s; c03r = c03;
	### Coeff d:
	d00 = - 2*c00r*x1 - c03r*b00r;
	d10 = 2*c00r - 2*c10r*x1 + 2*x1*b00r + s*b00r - c03r*b10r;
	d20 = 2*c10r - 2*c20r*x1 - 5*b00r + 2*x1*b10r + s*b10r - c03r*b20r;
	d30 = 2*c20r - c03r*x1 - 2*c30r*x1 + c03r*s - 5*b10r + 2*x1*b20r + s*b20r;
	d40 = 2*c30r + 2*x1^2 - x1*s - s^2 - 5*b20r;
	d50 = 5*(s - x1);
	# x3:
	d01 = - (2*c01r*x1 + c03r*b01r);
	d02 = - (2*c02r*x1 + c03r*b02r);
	d11 = 2*c01r - 2*c11r*x1 + 2*x1*b01r + s*b01r - c03r*b11r;
	d12 = 2*c02r - 2*c12r*x1 + 2*c03r*x1 + 3*c03r*s + 2*x1*b02r + s*b02r;
	d21 = 2*c11r - 2*c21r*x1 + c03r*x1 + 4*c03r*s - 5*b01r + 2*x1*b11r + s*b11r;
	d22 = 2*c12r - 4*c03r - 2*x1^2 + 4*x1*s - 3*s^2 - 5*b02r;
	d31 = 2*c21r - 3*c03r - 6*x1^2 - 5*x1*s - 4*s^2 - 5*b11r;
	#
	d01r = d01 - 15*x1*e4 - 4*s*e4;
	d11r = d11 + 15*e4 + 15*x1*e3 + 4*s*e3;
	d02r = d02 + 8*e4; d12r = d12 - 8*e3;
	d21r = d21 - 15*e3 - 15*x1*e2 - 4*s*e2;
	d31r = d31 + 15*e2 + 15*x1*s + 4*s^2;
	d22r = d22 + 8*e2;
	### Coeff x3:
	### x3: Div
	# x2^9
	dq9 = 48*x1^2 - 24*x1*s + 3*s^2 - 8*x1*d50 + 2*s*d50;
	dq8 = - 16*x1^3 - 56*x1^2*s + 31*x1*s^2 - 4*s^3 - 8*x1*d40 + 2*s*d40 + 8*x1^2*d50 - 2*x1*s*d50 +
		+ 24*x1*d22r - 6*s*d22r - 2*d50*d22r - 16*x1*d31r + 4*s*d31r;
	dq7 = 16*x1^2*b11r - 8*x1*s*b11r + s^2*b11r - 8*x1*d30 + 2*s*d30 + 8*x1^2*d40 - 2*x1*s*d40 +
		+ 24*x1*d12r - 6*s*d12r - 2*d50*d12r - 16*x1*d21r + 4*s*d21r - 8*x1^2*d22r - 30*x1*s*d22r +
		+ 8*s^2*d22r - 2*d40*d22r + 2*x1*d50*d22r + 3*d22r^2 + 8*x1^2*d31r + 10*x1*s*d31r - 3*s^2*d31r +
		- 4*d22r*d31r + 2*d31r^2;
	dq6 = 16*x1^2*b01r - 8*x1*s*b01r + s^2*b01r - 8*x1*d20 + 2*s*d20 + 8*x1^2*d30 - 2*x1*s*d30 +
		+ 24*x1*d02r - 6*s*d02r - 2*d50*d02r - 16*x1*d11r + 4*s*d11r - 8*x1^2*d12r - 30*x1*s*d12r +
		+ 8*s^2*d12r - 2*d40*d12r + 2*x1*d50*d12r + 8*x1^2*d21r + 10*x1*s*d21r - 3*s^2*d21r +
		+ 8*x1*b11r*d22r - 2*s*b11r*d22r - 2*d30*d22r + 2*x1*d40*d22r + 6*d12r*d22r - 4*d21r*d22r +
		- x1*d22r^2 - 4*s*d22r^2 - 4*x1*b02r*d31r + s*b02r*d31r - 4*d12r*d31r + 4*d21r*d31r +
		+ 2*x1*d22r*d31r + 3*s*d22r*d31r - 2*x1*d31r^2;
	dq5 = - 8*x1*d10 + 2*s*d10 + 8*x1^2*d20 - 2*x1*s*d20 - 16*x1*d01r + 4*s*d01r - 8*x1^2*d02r +
		- 30*x1*s*d02r + 8*s^2*d02r - 2*d40*d02r + 2*x1*d50*d02r + 8*x1^2*d11r + 10*x1*s*d11r +
		- 3*s^2*d11r + 8*x1*b11r*d12r - 2*s*b11r*d12r - 2*d30*d12r + 2*x1*d40*d12r + 3*d12r^2 +
		- 4*x1*b02r*d21r + s*b02r*d21r - 4*d12r*d21r + 2*d21r^2 + 8*x1*b01r*d22r - 2*s*b01r*d22r +
		- 2*d20*d22r + 2*x1*d30*d22r + 6*d02r*d22r - 4*d11r*d22r - 2*x1*d12r*d22r - 8*s*d12r*d22r +
		+ 2*x1*d21r*d22r + 3*s*d21r*d22r + b11r*d22r^2 - 4*d02r*d31r + 4*d11r*d31r + 2*x1*d12r*d31r +
		+ 3*s*d12r*d31r - 4*x1*d21r*d31r - b02r*d22r*d31r;
	dq4 = - 8*x1*d00 + 2*s*d00 + 8*x1^2*d10 - 2*x1*s*d10 + 8*x1^2*d01r + 10*x1*s*d01r - 3*s^2*d01r +
		+ 8*x1*b11r*d02r - 2*s*b11r*d02r - 2*d30*d02r + 2*x1*d40*d02r - 4*x1*b02r*d11r + s*b02r*d11r +
		+ 8*x1*b01r*d12r - 2*s*b01r*d12r - 2*d20*d12r + 2*x1*d30*d12r + 6*d02r*d12r - 4*d11r*d12r +
		- x1*d12r^2 - 4*s*d12r^2 - 4*d02r*d21r + 4*d11r*d21r + 2*x1*d12r*d21r + 3*s*d12r*d21r +
		- 2*x1*d21r^2 - 2*d10*d22r + 2*x1*d20*d22r - 4*d01r*d22r - 2*x1*d02r*d22r - 8*s*d02r*d22r +
		+ 2*x1*d11r*d22r + 3*s*d11r*d22r + 2*b11r*d12r*d22r - b02r*d21r*d22r + b01r*d22r^2 + 4*d01r*d31r +
		+ 2*x1*d02r*d31r + 3*s*d02r*d31r - 4*x1*d11r*d31r - b02r*d12r*d31r;
	dq3 = 8*x1^2*d00 - 2*x1*s*d00 - 4*x1*b02r*d01r + s*b02r*d01r + 8*x1*b01r*d02r - 2*s*b01r*d02r +
		- 2*d20*d02r + 2*x1*d30*d02r + 3*d02r^2 - 4*d02r*d11r + 2*d11r^2 - 2*d10*d12r + 2*x1*d20*d12r +
		- 4*d01r*d12r - 2*x1*d02r*d12r - 8*s*d02r*d12r + 2*x1*d11r*d12r + 3*s*d11r*d12r + b11r*d12r^2 +
		+ 4*d01r*d21r + 2*x1*d02r*d21r + 3*s*d02r*d21r - 4*x1*d11r*d21r - b02r*d12r*d21r - 2*d00*d22r +
		+ 2*x1*d10*d22r + 2*x1*d01r*d22r + 3*s*d01r*d22r + 2*b11r*d02r*d22r - b02r*d11r*d22r +
		+ 2*b01r*d12r*d22r - 4*x1*d01r*d31r - b02r*d02r*d31r;
	dq2 = - 2*d10*d02r + 2*x1*d20*d02r - 4*d01r*d02r - x1*d02r^2 - 4*s*d02r^2 + 4*d01r*d11r +
		+ 2*x1*d02r*d11r + 3*s*d02r*d11r - 2*x1*d11r^2 - 2*d00*d12r + 2*x1*d10*d12r + 2*x1*d01r*d12r +
		+ 3*s*d01r*d12r + 2*b11r*d02r*d12r - b02r*d11r*d12r + b01r*d12r^2 - 4*x1*d01r*d21r +
		- b02r*d02r*d21r + 2*x1*d00*d22r - b02r*d01r*d22r + 2*b01r*d02r*d22r;
	dq1 = 2*d01r^2 - 2*d00*d02r + 2*x1*d10*d02r + 2*x1*d01r*d02r + 3*s*d01r*d02r + b11r*d02r^2 +
		- 4*x1*d01r*d11r - b02r*d02r*d11r + 2*x1*d00*d12r - b02r*d01r*d12r + 2*b01r*d02r*d12r;
	dq0 = - 2*x1*d01r^2 + 2*x1*d00*d02r - b02r*d01r*d02r + b01r*d02r^2;
	
	### x3: x0
	# x2^10
	v10 = - 16*x1*d50 + 4*s*d50;
	v9 = 16*x1^3 - 24*x1^2*s + 9*x1*s^2 - s^3 - 16*x1*d40 + 4*s*d40 + 8*x1^2*d50 + 10*x1*s*d50 +
		- 3*s^2*d50 - 4*d50*d22r + 2*d50*d31r;
	v8 = 16*b20r*x1^2 - 8*b20r*x1*s + b20r*s^2 - 16*x1*d30 + 4*s*d30 + 8*x1^2*d40 + 10*x1*s*d40 +
		- 3*s^2*d40 - 4*x1*b02r*d50 + s*b02r*d50 - 4*d50*d12r + 2*d50*d21r + 8*x1^2*d22r - 10*x1*s*d22r +
		+ 2*s^2*d22r - 4*d40*d22r + 2*x1*d50*d22r + 3*s*d50*d22r + 2*d40*d31r - 2*x1*d50*d31r;
	v7 = 16*b10r*x1^2 - 8*b10r*x1*s + b10r*s^2 - 16*x1*d20 + 4*s*d20 + 8*x1^2*d30 + 10*x1*s*d30 +
		- 3*s^2*d30 - 4*x1*b02r*d40 + s*b02r*d40 - 4*d50*d02r + 2*d50*d11r + 8*x1^2*d12r - 10*x1*s*d12r +
		+ 2*s^2*d12r - 4*d40*d12r + 2*x1*d50*d12r + 3*s*d50*d12r + 2*d40*d21r - 2*x1*d50*d21r +
		+ 8*b20r*x1*d22r - 2*b20r*s*d22r - 4*d30*d22r + 2*x1*d40*d22r + 3*s*d40*d22r - b02r*d50*d22r +
		+ x1*d22r^2 - s*d22r^2 + 2*d30*d31r - 2*x1*d40*d31r;
	v6 = 16*b00r*x1^2 - 8*b00r*x1*s + b00r*s^2 - 16*x1*d10 + 4*s*d10 + 8*x1^2*d20 + 10*x1*s*d20 +
		- 3*s^2*d20 - 4*x1*b02r*d30 + s*b02r*d30 + 2*d50*d01r + 8*x1^2*d02r - 10*x1*s*d02r + 2*s^2*d02r +
		- 4*d40*d02r + 2*x1*d50*d02r + 3*s*d50*d02r + 2*d40*d11r - 2*x1*d50*d11r + 8*b20r*x1*d12r +
		- 2*b20r*s*d12r - 4*d30*d12r + 2*x1*d40*d12r + 3*s*d40*d12r - b02r*d50*d12r + 2*d30*d21r +
		- 2*x1*d40*d21r + 8*b10r*x1*d22r - 2*b10r*s*d22r - 4*d20*d22r + 2*x1*d30*d22r + 3*s*d30*d22r +
		- b02r*d40*d22r + 2*x1*d12r*d22r - 2*s*d12r*d22r + b20r*d22r^2 + 2*d20*d31r - 2*x1*d30*d31r;
	v5 = - 16*x1*d00 + 4*s*d00 + 8*x1^2*d10 + 10*x1*s*d10 - 3*s^2*d10 - 4*x1*b02r*d20 + s*b02r*d20 +
		+ 2*d40*d01r - 2*x1*d50*d01r + 8*b20r*x1*d02r - 2*b20r*s*d02r - 4*d30*d02r + 2*x1*d40*d02r +
		+ 3*s*d40*d02r - b02r*d50*d02r + 2*d30*d11r - 2*x1*d40*d11r + 8*b10r*x1*d12r - 2*b10r*s*d12r +
		- 4*d20*d12r + 2*x1*d30*d12r + 3*s*d30*d12r - b02r*d40*d12r + x1*d12r^2 - s*d12r^2 +
		+ 2*d20*d21r - 2*x1*d30*d21r + 8*b00r*x1*d22r - 2*b00r*s*d22r - 4*d10*d22r + 2*x1*d20*d22r +
		+ 3*s*d20*d22r - b02r*d30*d22r + 2*x1*d02r*d22r - 2*s*d02r*d22r + 2*b20r*d12r*d22r + b10r*d22r^2 +
		+ 2*d10*d31r - 2*x1*d20*d31r;
	v4 = 8*x1^2*d00 + 10*x1*s*d00 - 3*s^2*d00 - 4*x1*b02r*d10 + s*b02r*d10 + 2*d30*d01r - 2*x1*d40*d01r +
		+ 8*b10r*x1*d02r - 2*b10r*s*d02r - 4*d20*d02r + 2*x1*d30*d02r + 3*s*d30*d02r - b02r*d40*d02r +
		+ 2*d20*d11r - 2*x1*d30*d11r + 8*b00r*x1*d12r - 2*b00r*s*d12r - 4*d10*d12r + 2*x1*d20*d12r +
		+ 3*s*d20*d12r - b02r*d30*d12r + 2*x1*d02r*d12r - 2*s*d02r*d12r + b20r*d12r^2 + 2*d10*d21r +
		- 2*x1*d20*d21r - 4*d00*d22r + 2*x1*d10*d22r + 3*s*d10*d22r - b02r*d20*d22r + 2*b20r*d02r*d22r +
		+ 2*b10r*d12r*d22r + b00r*d22r^2 + 2*d00*d31r - 2*x1*d10*d31r;
	v3 = - 4*x1*b02r*d00 + s*b02r*d00 + 2*d20*d01r - 2*x1*d30*d01r + 8*b00r*x1*d02r - 2*b00r*s*d02r +
		- 4*d10*d02r + 2*x1*d20*d02r + 3*s*d20*d02r - b02r*d30*d02r + x1*d02r^2 - s*d02r^2 + 2*d10*d11r +
		- 2*x1*d20*d11r - 4*d00*d12r + 2*x1*d10*d12r + 3*s*d10*d12r - b02r*d20*d12r + 2*b20r*d02r*d12r +
		+ b10r*d12r^2 + 2*d00*d21r - 2*x1*d10*d21r + 2*x1*d00*d22r + 3*s*d00*d22r - b02r*d10*d22r +
		+ 2*b10r*d02r*d22r + 2*b00r*d12r*d22r - 2*x1*d00*d31r;
	v2 = 2*d10*d01r - 2*x1*d20*d01r - 4*d00*d02r + 2*x1*d10*d02r + 3*s*d10*d02r - b02r*d20*d02r +
		+ b20r*d02r^2 + 2*d00*d11r - 2*x1*d10*d11r + 2*x1*d00*d12r + 3*s*d00*d12r - b02r*d10*d12r +
		+ 2*b10r*d02r*d12r + b00r*d12r^2 - 2*x1*d00*d21r - b02r*d00*d22r + 2*b00r*d02r*d22r;
	v1 = 2*d00*d01r - 2*x1*d10*d01r + 2*x1*d00*d02r + 3*s*d00*d02r - b02r*d10*d02r + b10r*d02r^2 +
		- 2*x1*d00*d11r - b02r*d00*d12r + 2*b00r*d02r*d12r;
	v0 = - 2*x1*d00*d01r - b02r*d00*d02r + b00r*d02r^2;
	# Reduced
	# x2^3
	dq3r = dq3 - dq7*e4 + dq6*e3 + dq9*e3^2 - dq5*e2 + 2*dq9*e4*e2 - 2*dq8*e3*e2 + dq7*e2^2 +
		- dq9*e2^3 + dq4*s - 2*dq8*e4*s + 2*dq7*e3*s - 2*dq6*e2*s - 6*dq9*e3*e2*s + 3*dq8*e2^2*s +
		+ dq5*s^2 - 3*dq9*e4*s^2 + 3*dq8*e3*s^2 - 3*dq7*e2*s^2 + 6*dq9*e2^2*s^2 + dq6*s^3 +
		+ 4*dq9*e3*s^3 - 4*dq8*e2*s^3 + dq7*s^4 - 5*dq9*e2*s^4 + dq8*s^5 + dq9*s^6;
	dq2r = dq2 - dq6*e4 + dq5*e3 - 2*dq9*e4*e3 + dq8*e3^2 - dq4*e2 + 2*dq8*e4*e2 - 2*dq7*e3*e2 +
		+ dq6*e2^2 + 3*dq9*e3*e2^2 - dq8*e2^3 - dq7*e4*s + dq6*e3*s + 2*dq9*e3^2*s - dq5*e2*s +
		+ 4*dq9*e4*e2*s - 4*dq8*e3*e2*s + 2*dq7*e2^2*s - 3*dq9*e2^3*s - dq8*e4*s^2 + dq7*e3*s^2 +
		- dq6*e2*s^2 - 6*dq9*e3*e2*s^2 + 3*dq8*e2^2*s^2 - dq9*e4*s^3 + dq8*e3*s^3 - dq7*e2*s^3 +
		+ 4*dq9*e2^2*s^3 + dq9*e3*s^4 - dq8*e2*s^4 - dq9*e2*s^5;
	dq1r = dq1 - dq5*e4 + dq9*e4^2 + dq4*e3 - 2*dq8*e4*e3 + dq7*e3^2 + dq7*e4*e2 - dq6*e3*e2 +
		- 2*dq9*e3^2*e2 - dq9*e4*e2^2 + dq8*e3*e2^2 - dq6*e4*s + dq5*e3*s - 4*dq9*e4*e3*s +
		+ 2*dq8*e3^2*s + 2*dq8*e4*e2*s - 2*dq7*e3*e2*s + 3*dq9*e3*e2^2*s - dq7*e4*s^2 + dq6*e3*s^2 +
		+ 3*dq9*e3^2*s^2 + 3*dq9*e4*e2*s^2 - 3*dq8*e3*e2*s^2 - dq8*e4*s^3 + dq7*e3*s^3 +
		- 4*dq9*e3*e2*s^3 - dq9*e4*s^4 + dq8*e3*s^4 + dq9*e3*s^5;
	dq0r = dq0 - dq4*e4 + dq8*e4^2 - dq7*e4*e3 + dq6*e4*e2 + 2*dq9*e4*e3*e2 - dq8*e4*e2^2 +
		- dq5*e4*s + 2*dq9*e4^2*s - 2*dq8*e4*e3*s + 2*dq7*e4*e2*s - 3*dq9*e4*e2^2*s - dq6*e4*s^2 +
		- 3*dq9*e4*e3*s^2 + 3*dq8*e4*e2*s^2 - dq7*e4*s^3 + 4*dq9*e4*e2*s^3 - dq8*e4*s^4 - dq9*e4*s^5;
	#
	v3r = v3 - v7*e4 + v6*e3 - 2*v10*e4*e3 + v9*e3^2 - v5*e2 + 2*v9*e4*e2 - 2*v8*e3*e2 + v7*e2^2 +
		+ 3*v10*e3*e2^2 - v9*e2^3 + v4*s - 2*v8*e4*s + 2*v7*e3*s + 3*v10*e3^2*s - 2*v6*e2*s +
		+ 6*v10*e4*e2*s - 6*v9*e3*e2*s + 3*v8*e2^2*s - 4*v10*e2^3*s + v5*s^2 - 3*v9*e4*s^2 +
		+ 3*v8*e3*s^2 - 3*v7*e2*s^2 - 12*v10*e3*e2*s^2 + 6*v9*e2^2*s^2 + v6*s^3 - 4*v10*e4*s^3 +
		+ 4*v9*e3*s^3 - 4*v8*e2*s^3 + 10*v10*e2^2*s^3 + v7*s^4 + 5*v10*e3*s^4 - 5*v9*e2*s^4 +
		+ v8*s^5 - 6*v10*e2*s^5 + v9*s^6 + v10*s^7;
	v2r = v2 - v6*e4 + v10*e4^2 + v5*e3 - 2*v9*e4*e3 + v8*e3^2 - v4*e2 + 2*v8*e4*e2 - 2*v7*e3*e2 +
		- 3*v10*e3^2*e2 + v6*e2^2 - 3*v10*e4*e2^2 + 3*v9*e3*e2^2 - v8*e2^3 + v10*e2^4 - v7*e4*s +
		+ v6*e3*s - 4*v10*e4*e3*s + 2*v9*e3^2*s - v5*e2*s + 4*v9*e4*e2*s - 4*v8*e3*e2*s + 2*v7*e2^2*s +
		+ 9*v10*e3*e2^2*s - 3*v9*e2^3*s - v8*e4*s^2 + v7*e3*s^2 + 3*v10*e3^2*s^2 - v6*e2*s^2 +
		+ 6*v10*e4*e2*s^2 - 6*v9*e3*e2*s^2 + 3*v8*e2^2*s^2 - 6*v10*e2^3*s^2 - v9*e4*s^3 + v8*e3*s^3 +
		- v7*e2*s^3 - 8*v10*e3*e2*s^3 + 4*v9*e2^2*s^3 - v10*e4*s^4 + v9*e3*s^4 - v8*e2*s^4 +
		+ 5*v10*e2^2*s^4 + v10*e3*s^5 - v9*e2*s^5 - v10*e2*s^6;
	v1r = v1 - v5*e4 + v9*e4^2 + v4*e3 - 2*v8*e4*e3 + v7*e3^2 + v10*e3^3 + v7*e4*e2 - v6*e3*e2 +
		+ 4*v10*e4*e3*e2 - 2*v9*e3^2*e2 - v9*e4*e2^2 + v8*e3*e2^2 - v10*e3*e2^3 - v6*e4*s +
		+ 2*v10*e4^2*s + v5*e3*s - 4*v9*e4*e3*s + 2*v8*e3^2*s + 2*v8*e4*e2*s - 2*v7*e3*e2*s +
		- 6*v10*e3^2*e2*s - 3*v10*e4*e2^2*s + 3*v9*e3*e2^2*s - v7*e4*s^2 + v6*e3*s^2 - 6*v10*e4*e3*s^2 +
		+ 3*v9*e3^2*s^2 + 3*v9*e4*e2*s^2 - 3*v8*e3*e2*s^2 + 6*v10*e3*e2^2*s^2 - v8*e4*s^3 +
		+ v7*e3*s^3 + 4*v10*e3^2*s^3 + 4*v10*e4*e2*s^3 - 4*v9*e3*e2*s^3 - v9*e4*s^4 + v8*e3*s^4 +
		- 5*v10*e3*e2*s^4 - v10*e4*s^5 + v9*e3*s^5 + v10*e3*s^6;
	v0r = v0 - v4*e4 + v8*e4^2 - v7*e4*e3 - v10*e4*e3^2 + v6*e4*e2 - 2*v10*e4^2*e2 + 2*v9*e4*e3*e2 +
		- v8*e4*e2^2 + v10*e4*e2^3 - v5*e4*s + 2*v9*e4^2*s - 2*v8*e4*e3*s + 2*v7*e4*e2*s +
		+ 6*v10*e4*e3*e2*s - 3*v9*e4*e2^2*s - v6*e4*s^2 + 3*v10*e4^2*s^2 - 3*v9*e4*e3*s^2 +
		+ 3*v8*e4*e2*s^2 - 6*v10*e4*e2^2*s^2 - v7*e4*s^3 - 4*v10*e4*e3*s^3 + 4*v9*e4*e2*s^3 +
		- v8*e4*s^4 + 5*v10*e4*e2*s^4 - v9*e4*s^5 - v10*e4*s^6;
	###
	coeff = list(
		v  = cbind(v0r, v1r, v2r, v3r),
		dq = cbind(dq0r, dq1r, dq2r, dq3r));
	return(coeff);
}

