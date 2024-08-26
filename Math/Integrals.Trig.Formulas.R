
### Trig: Formulas

### Trig(3*x)
x = pi/7
cos(3*x) - 4*cos(x)^3 + 3*cos(x)
sin(3*x) + 4*sin(x)^3 - 3*sin(x)
#
tan(3*x) + (4*sin(x)^3 - 3*sin(x)) / (4*cos(x)^3 - 3*cos(x))
tan(3*x) + (tan(x)^3 - 3*tan(x)) / (1 - 3*tan(x)^2)

### ATAN
x = sqrt(5)
# Note: force result between [-pi/2, pi/2];
atan((3*x - x^3) / (1 - 3*x^2))
3 * atan(x) - pi
#
atan((1 - 3*x^2) / (3*x - x^3))
3/2*pi - 3 * atan(x)

# Transforms:
atan((2*x^3 - 3*x + 1) / (2*x^3 - 6*x^2 + 3*x))
- 3 * atan(1 - 1/x)
# TODO: more;


### 1/3 Angle:
x = pi/7
cos(x/3) # ==
((cos(x) + 1i*sin(x))^(1/3) + (cos(x) - 1i*sin(x))^(1/3)) / 2

### Note:
# D(cos(3*x)) => - 3*sin(3*x), but masked;
x = pi/7
- 12*sin(x)*cos(x)^2 + 3*sin(x)
- 3*sin(3*x)


### Examples:
x = pi/7
integrate(\(x) -(4*sin(x)^3 - 3*sin(x)) / (4*cos(x)^3 - 3*cos(x)), 0, up=x)
integrate(\(x) -(tan(x)^3 - 3*tan(x)) / (1 - 3*tan(x)^2), 0, up=x)
integrate(\(x) -(1 - 4*x^2) / (4*x^3 - 3*x), lower = cos(x), 1)
integrate(\(x) -(x^3 - 3*x) / ((x^2 + 1) * (1 - 3*x^2)), 0, up = tan(x))
integrate(\(x) tan(3*x), 0, up=x)
- log(cos(3*x)) / 3;

