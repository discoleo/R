###################
##
## Mixed Fractions: Trig & Polynomial
## Extended Variants


### I( x / (b - sin(x)) )
### I( x / (b + sin(x)) )
# Maths 505: Thank you for this wonderful integral
# https://www.youtube.com/watch?v=zaiA-0OBKa8
# Maths 505: SUPREME GOLDEN INTEGRAL
# https://www.youtube.com/watch?v=uXHBzyypGOw

b = sqrt(5) # b > 1;
integrate(\(x) x / (b - sin(x)), 0, pi)
pi / sqrt(b^2 - 1) * (pi - atan(sqrt(b^2 - 1)))

integrate(\(x) x / (b + sin(x)), 0, pi)
pi / sqrt(b^2 - 1) * atan(sqrt(b^2 - 1))

integrate(\(x) x / (b - sin(x)), 0, 2*pi)
pi / sqrt(b^2 - 1) * (pi + 2*atan(sqrt(b^2 - 1)))

integrate(\(x) x / (b + sin(x)), 0, 2*pi)
pi / sqrt(b^2 - 1) * (3*pi - 2*atan(sqrt(b^2 - 1)))

### Helper
integrate(\(x) 1 / (b - sin(x)), 0, pi)
integrate(\(x) 1 / (b + sin(x)), pi, 2*pi)
2 / sqrt(b^2 - 1) * (pi - atan(sqrt(b^2 - 1)))

integrate(\(x) 1 / (b + sin(x)), 0, pi)
integrate(\(x) 1 / (b - sin(x)), pi, 2*pi)
2 / sqrt(b^2 - 1) * atan(sqrt(b^2 - 1))


### More Helper:
b = sqrt(5)
integrate(\(x) x / (b - sin(x)), pi, 2*pi)
3*pi / sqrt(b^2 - 1) * atan(sqrt(b^2 - 1))

integrate(\(x) x / (b + sin(x)), pi, 2*pi)
3*pi / sqrt(b^2 - 1) * (pi - atan(sqrt(b^2 - 1)))


### Gen 1:
b = sqrt(5)
integrate(\(x) x / (b - sin(x)), 2*pi, 3*pi)
5*pi / sqrt(b^2 - 1) * (pi - atan(sqrt(b^2 - 1)))

integrate(\(x) x / (b + sin(x)), 2*pi, 3*pi)
5*pi / sqrt(b^2 - 1) * atan(sqrt(b^2 - 1))


### Gen 2: cos-Variants
b = sqrt(7)
integrate(\(x) 1 / (b - cos(x)), 0, pi)
integrate(\(x) 1 / (b + cos(x)), 0, pi)
pi / sqrt(b^2 - 1)

#
b = sqrt(7)
integrate(\(x) 1 / (b^2 - cos(x)^2), 0, pi)
pi / (b * sqrt(b^2 - 1))
#
integrate(\(x) 1 / (b^2 + cos(x)^2), 0, pi)
pi / (b * sqrt(b^2 + 1))
# Note: Power can be extended;


# I( x / (b^2 - cos(x)^2) )
b = sqrt(5)
integrate(\(x) x / (b^2 - cos(x)^2), 0, pi)
pi^2 / (2*b*sqrt(b^2 - 1))

# I( x / (b^2 + cos(x)^2) )
integrate(\(x) x / (b^2 + cos(x)^2), 0, pi)
pi^2 / (2*b*sqrt(b^2 + 1))

# on [pi, 2*pi]
integrate(\(x) x / (b^2 - cos(x)^2), pi, 2*pi)
3*pi^2 / (2*b * sqrt(b^2 - 1))
#
integrate(\(x) x / (b^2 + cos(x)^2), pi, 2*pi)
3*pi^2 / (2*b * sqrt(b^2 + 1))


# Cases:
integrate(\(x) x / (4 - cos(x)^2), 0, pi)
pi^2 / (4*sqrt(3))

integrate(\(x) x / (9 - cos(x)^2), 0, pi)
pi^2 / (6*sqrt(8))


###
b = sqrt(5);
integrate(\(x) x * cos(x) / (b^2 - cos(x)^2), 0, pi)
- 1/sqrt(b^2 - 1) * integrate(\(x) atan(sin(x) / sqrt(b^2 - 1)), 0, pi)$value

# TODO


### Simple Variants:
b = sqrt(5); up = pi / 5
integrate(\(x) cos(x) / (b^2 - cos(x)^2), 0, up)
integrate(\(x) 1 / (x^2 + b^2 - 1), 0, sin(up))
1/sqrt(b^2 - 1) * atan(sin(up) / sqrt(b^2 - 1))
#
integrate(\(x) cos(x) / (b^2 + cos(x)^2), 0, up)
1/(2*sqrt(b^2 + 1)) * log((sqrt(b^2 + 1) + sin(up)) / (sqrt(b^2 + 1) - sin(up)))


###
b = sqrt(7)
integrate(\(x) atan(sin(x) / sqrt(b^2 - 1)), 0, pi)
# Weierstrass:
integrate(\(x) 2 * atan(2*x / (sqrt(b^2 - 1) * (x^2+1))) / (x^2 + 1), 0, Inf)

# Alternative:
# - but may be circular;
k = sqrt(3); # k = sqrt(3) + 1E-4
integrate(\(x) atan(k * sin(x) / sqrt(b^2 - 1)), 0, pi)
# dk =>
sqrt(b^2 - 1) * integrate(\(x) sin(x) / (k^2*sin(x)^2 + b^2 - 1), 0, pi)$value

# TODO
