


### Trig: TAN


### Helper Constants
Catalan = 0.915965594177219015054603514;


##################
##################

###
integrate(\(x) x * tan(x), 0, pi/4)
Catalan / 2 - 1/8 * pi*log(2)

### [by Parts]
integrate(\(x) x^2 * tan(x)^2, 0, pi/4)
- 1/3*(pi/4)^3 + pi^2 / 16 + 1/4 * pi*log(2) - Catalan


###
integrate(\(x) x^2 * tan(x), 0, pi/4)
- 21/64 * zeta(3) - 1/32 * pi^2 * log(2) + 1/4 * pi * Catalan

