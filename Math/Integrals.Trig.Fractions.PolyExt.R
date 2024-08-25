###################
##
## Mixed Fractions: Trig & Polynomial
## Extended Variants


### I( x / (b - sin(x)) )
### I( x / (b + sin(x)) )
# Maths 505: Thank you for this wonderful integral
# https://www.youtube.com/watch?v=zaiA-0OBKa8

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


### Gen 2:
b = sqrt(7)
integrate(\(x) 1 / (b - cos(x)), 0, pi)
integrate(\(x) 1 / (b + cos(x)), 0, pi)
pi / sqrt(b^2 - 1)

# TODO


#
b = sqrt(5)
integrate(\(x) x / (b^2 - cos(x)^2), 0, pi)
pi^2 / (2*b*sqrt(b^2 - 1))

integrate(\(x) x / (4 - cos(x)^2), 0, pi)
pi^2 / (4*sqrt(3))

integrate(\(x) x / (9 - cos(x)^2), 0, pi)
pi^2 / (6*sqrt(8))

