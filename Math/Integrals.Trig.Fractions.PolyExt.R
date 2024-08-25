###################
##
## Mixed Fractions: Trig & Polynomial
## Extended Variants


### I( x / (k - sin(x)) )
### I( x / (k + sin(x)) )
# Maths 505: Thank you for this wonderful integral
# https://www.youtube.com/watch?v=zaiA-0OBKa8

k = sqrt(5) # k > 1;
integrate(\(x) x / (k - sin(x)), 0, pi)
pi / sqrt(k^2 - 1) * (pi - atan(sqrt(k^2 - 1)))

integrate(\(x) x / (k + sin(x)), 0, pi)
pi / sqrt(k^2 - 1) * atan(sqrt(k^2 - 1))

integrate(\(x) x / (k - sin(x)), 0, 2*pi)
pi / sqrt(k^2 - 1) * (pi + 2*atan(sqrt(k^2 - 1)))

integrate(\(x) x / (k + sin(x)), 0, 2*pi)
pi / sqrt(k^2 - 1) * (3*pi - 2*atan(sqrt(k^2 - 1)))

### Helper
integrate(\(x) 1 / (k - sin(x)), 0, pi)
integrate(\(x) 1 / (k + sin(x)), pi, 2*pi)
2 / sqrt(k^2 - 1) * (pi - atan(sqrt(k^2 - 1)))

integrate(\(x) 1 / (k + sin(x)), 0, pi)
integrate(\(x) 1 / (k - sin(x)), pi, 2*pi)
2 / sqrt(k^2 - 1) * atan(sqrt(k^2 - 1))


### More Helper:
k = sqrt(5)
integrate(\(x) x / (k - sin(x)), pi, 2*pi)
3*pi / sqrt(k^2 - 1) * atan(sqrt(k^2 - 1))

integrate(\(x) x / (k + sin(x)), pi, 2*pi)
3*pi / sqrt(k^2 - 1) * (pi - atan(sqrt(k^2 - 1)))


### Gen 1:
k = sqrt(5)
integrate(\(x) x / (k - sin(x)), 2*pi, 3*pi)
5*pi / sqrt(k^2 - 1) * (pi - atan(sqrt(k^2 - 1)))

integrate(\(x) x / (k + sin(x)), 2*pi, 3*pi)
5*pi / sqrt(k^2 - 1) * atan(sqrt(k^2 - 1))

