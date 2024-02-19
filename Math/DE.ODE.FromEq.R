

### ODE From another Equation


################

### Basic Examples
# see file: DE.ODE.Polynomial.R


### 2*atan((y + 2) / (1 - x)) - log(y + 2) = Const
# Example from:
# Maths 505: A very interesting differential equation
# https://www.youtube.com/watch?v=zDw4FjmLv0Y

dy - 2*(y + 2)^2 / (y + x + 1)^2 = 0


##################

### y * exp(x*y) + x * exp(y^2) = F(x)

# for testing:
x = sqrt(5); y = x; f = 2*x*exp(x^2); dy = 1; d2y = 0;
df = 2*(1 + 2*x^2)*exp(x^2); d2f = 2*(4*x + 2*x*(1 + 2*x^2))*exp(x^2);

# D =>
(dy + x*y*dy + y^2) * exp(x*y) + (1 + 2*x*y*dy) * exp(y^2) - df # = 0
# x*exp(y^2) = - y * exp(x*y) + f;
x*(dy + x*y*dy + y^2) * exp(x*y) - (2*x*y*dy + 1) * (y * exp(x*y) - f) - x*df # = 0
(x*dy - 2*x*y^2*dy + x^2*y*dy + x*y^2 - y) * exp(x*y) + 2*x*f*y*dy + f - x*df # = 0

# D2 =>
(x^2*y*d2y - 2*x*y^2*d2y + x*d2y - 2*x^2*y^2*dy^2 + x^3*y*dy^2 - 4*x*y*dy^2 + 2*x^2*dy^2 +
	- 2*x*y^3*dy + 2*x^2*y^2*dy - 2*y^2*dy + 2*x*y*dy + 2*x*y*dy + x*y^3) * exp(x*y) +
	+ 2*x*f*y*d2y + 2*x*f*dy^2 + 2*f*y*dy + 2*x*df*y*dy - x*d2f # = 0


# TODO:
# - substitute exp(x*y) and check;

