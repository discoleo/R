

### Dogbone Integration

# 1. Complex Analysis: Dogbone Contour Example
# https://www.youtube.com/watch?v=UDIKojCQ94U
integrate(function(x) x^(3/4) * (3 - x)^(1/4) / (5 - x), lower=0, upper=3)

# 2. Complex Analysis: Dogbone Contour Example #2
# https://www.youtube.com/watch?v=q1BxM1MWAqA
integrate(function(x) (1 - x)^(-2/3) * (1 + x)^(-1/3) / (4 + x^2), lower=-1, upper=1)


# 3. Complex Analysis: Dogbone Contour Example #3
# https://www.youtube.com/watch?v=-HWcFun7e4k
# - is similar to [2];
integrate(function(x) (1 - x)^(1/2) * (1 + x)^(1/2) / (1 + x^2), lower=-1, upper=1)


# 4. qncubed3: Complex Analysis: Dogbone Contour Generalisation
# https://www.youtube.com/watch?v=w-NIlyXZzqU
p = 1/3; a = 1; b = 4;
integrate(function(x) (x - a)^p * (b - x)^(1 - p), lower=a, upper=b)

# 5. qncubed3: Complex Analysis: Double Keyhole Contour
# https://www.youtube.com/watch?v=LT3jpvWMH2s
integrate(function(x) 1 / (x*sqrt(x^2 - 1)), lower=1, upper=Inf)



##############

### Example 1:
# Complex Analysis: Dogbone Contour Example
# https://www.youtube.com/watch?v=UDIKojCQ94U

integrate(function(x) x^(3/4)*(3 - x)^(1/4) / (5 - x), lower=0, upper=3)
pi*(17/4 - (5^3*2)^(1/4))*sqrt(2)


### Generalizations:

integrate(function(x) x^(3/4)*(3 - x)^(1/4) / (7 - x), lower=0, upper=3)
pi*(7 - 3/4 - (7^3*4)^(1/4))*sqrt(2)

### Gen 1:
k = 11
integrate(function(x) x^(3/4)*(3 - x)^(1/4) / (k - x), lower=0, upper=3)
pi*(k - 3/4 - (k^3*(k - 3))^(1/4))*sqrt(2)

### Gen 2:
k = 11; b = 4;
integrate(function(x) x^(3/4)*(b - x)^(1/4) / (k - x), lower=0, upper=b)
pi*(k - b/4 - (k^3*(k - b))^(1/4))*sqrt(2)

### Gen 3: Full
k = 11; b = 4;
p = 5;
integrate(function(x) x^((p-1)/p)*(b - x)^(1/p) / (k - x), lower=0, upper=b)
2i * pi*(k - b/p - (k^(p-1)*(k - b))^(1/p)) * exp(1i*pi/p) / (exp(2i*pi/p) - 1)
pi*(k - b/p - (k^(p-1)*(k - b))^(1/p)) / sin(pi/p)

# Arg - Inf from above: (p+1)/p * pi;
# Arg - Inf from below: -(p-1)/p * pi;
# Continuous: OK
# 2*pi - (p-1)/p * pi => (p+1)/p * pi;

# Note:
# - variation/partitioning of powers can be emulated with fractional p;

### Ex:
p = 5/2;
k = 11; b = 4;
integrate(function(x) x^(3/5)*(b - x)^(2/5) / (k - x), lower=0, upper=b)
#
integrate(function(x) x^((p-1)/p)*(b - x)^(1/p) / (k - x), lower=0, upper=b)
pi*(k - b/p - (k^(p-1)*(k - b))^(1/p)) / sin(pi/p)


######################

### Example 2:
# Complex Analysis: Dogbone Contour Example #2
# https://www.youtube.com/watch?v=q1BxM1MWAqA

integrate(function(x) (1 - x)^(-2/3)*(1 + x)^(-1/3) / (4 + x^2), lower=-1, upper=1)
pi*sin(atan(2)/3 + pi/3) / sin(pi/3) / sqrt(5) / 2


### Gen 1:
k = sqrt(5)
integrate(function(x) (1 - x)^(-2/3)*(1 + x)^(-1/3) / (k^2 + x^2), lower=-1, upper=1)
pi*sin(atan(k)/3 + pi/3) / sin(pi/3) / sqrt(k^2 + 1) / k


### Gen 2:
k = sqrt(3);
p = 5;
integrate(function(x) (1 - x)^(1/p - 1) * (1 + x)^(-1/p) / (k^2 + x^2), lower=-1, upper=1)
pi*sin(atan(k)*(1-2/p) + pi/p) / sin(pi/p) / sqrt(1 + k^2) / k


### Gen 3:
b = 5;
k = sqrt(3);
p = 5;
integrate(function(x) (b - x)^(1/p - 1) * (b + x)^(-1/p) / (k^2 + x^2), lower=-b, upper=b)
pi*sin(atan(k/b)*(1-2/p) + pi/p) / sin(pi/p) / sqrt(b^2 + k^2) / k


### Derivation:
pi*(exp(1i*atan(k)*(1-2/p) + 2i*pi/p) +
	- exp(1i*(2/p-1)*atan(k))) / sqrt(1 + k^2) / (2i*k) / (exp(1i*pi/p) * sin(pi/p))

# Arg - Inf from above: (1/p - 2)*pi;
# Arg - Inf from below: 1/p * pi;
# Continuous: OK
# 2*pi + (1/p - 2)*pi => 1/p * pi;

# (1/p - 1)*2*pi (equivalent to: 2*pi/p); 0;
# (2*pi - fi)*(1/p-1) - fi/p = 2*pi/p + (1 - 2/p)*fi;


### Variant: n = 0
# (b - x)^(1/p - 0)
# Conditions: b > 0
b = sqrt(3);
k = sqrt(7);
p = 5;
integrate(function(x) (b - x)^(1/p) * (b + x)^(-1/p) / (k^2 + x^2), lower=-b, upper=b)
pi*sin(atan(k/b)*(-2/p) + pi/p) / sin(pi/p) / k


### Variant: n1 = -1; n2 = 1;
# - equivalent to substitution:
#   1/p => 1/p - 1 (in Variant 1);
b = sqrt(3);
k = sqrt(7);
p = 5;
integrate(function(x) (b - x)^(1/p - 1) * (b + x)^(1-1/p) / (k^2 + x^2), lower=-b, upper=b)
pi*sin(atan(k/b)*(2-2/p) + pi/p) / sin(pi/p) / k


#########################
#########################

# 4. qncubed3: Complex Analysis: Dogbone Contour Generalisation
# https://www.youtube.com/watch?v=w-NIlyXZzqU

p = 1/3; a = 1; b = 4;
integrate(function(x) (x - a)^p * (b - x)^(1 - p), lower=a, upper=b)
pi/2 * p*(1 - p)*(b - a)^2/sin(pi*p)


### Gen 1:
p = sqrt(3) - sqrt(2)
a = sqrt(5); b = sqrt(7);
integrate(function(x) x * (x - a)^p * (b - x)^(1 - p), lower=a, upper=b)
pi/6 * p*(1-p)*((2-p)*a + (1+p)*b)*(b - a)^2/sin(pi*p)


### Gen 2 & 3: (x - a)^p * (b - x)^q
# Transform: y = b - (b-a)/(x+1);
# New: (b-a)^(p+q+1) * x^p / (x+1)^(p+q+2);
# Interval: [0, 1] => [0, Inf];

### Gen 2:
p = 1/3; a = 1; b = 4;
integrate(function(x) (x - a)^p * (b - x)^p, lower=a, upper=b)
(b-a)^(2*p+1) * gamma(p+1)*gamma(p+1) / gamma(2*p+2)

###
p = sqrt(5) - sqrt(3); a = sqrt(2); b = 4;
integrate(function(x) (x - a)^p * (b - x)^p, lower=a, upper=b)
(b-a)^(2*p+1) * gamma(p+1)*gamma(p+1) / gamma(2*p+2)


### Gen 3:
p = 1/3; q = 1/5; a = 1; b = 4;
integrate(function(x) (x - a)^p * (b - x)^q, lower=a, upper=b)
(b-a)^(p+q+1) * gamma(p+1)*gamma(q+1) / gamma(p+q+2)

###
p = sqrt(3); q = sqrt(5); a = 1; b = 4;
integrate(function(x) (x - a)^p * (b - x)^q, lower=a, upper=b)
(b-a)^(p+q+1) * gamma(p+1)*gamma(q+1) / gamma(p+q+2)


### Variant:
p = sqrt(3); q = sqrt(5); a = 1; b = 4;
integrate(function(x) x * (x - a)^p * (b - x)^q, lower=a, upper=b)
b*(b-a)^(p+q+1) * gamma(p+1)*gamma(q+1) / gamma(p+q+2) +
	- (b-a)^(p+q+2) * gamma(p+1)*gamma(q+2) / gamma(p+q+3)
# simplified:
(b*(p+q+2) - (b-a)*(q+1)) * (b-a)^(p+q+1) * gamma(p+1)*gamma(q+1) / gamma(p+q+3)
(b*(p+1) + a*(q+1)) * (b-a)^(p+q+1) * gamma(p+1)*gamma(q+1) / gamma(p+q+3)


### Variant:
p = sqrt(3); q = sqrt(5); a = 1; b = 4;
integrate(function(x) x^2 * (x - a)^p * (b - x)^q, lower=a, upper=b)
b^2*(b-a)^(p+q+1) * gamma(p+1)*gamma(q+1) / gamma(p+q+2) +
	- 2*b*(b-a)^(p+q+2) * gamma(p+1)*gamma(q+2) / gamma(p+q+3) +
	+ (b-a)^(p+q+3) * gamma(p+1)*gamma(q+3) / gamma(p+q+4)
# simplified:
(b^2*(p+q+2)*(p+q+3) - 2*b*(b-a)*(q+1)*(p+q+3) + (b-a)^2*(q+1)*(q+2)) *
	(b-a)^(p+q+1) * gamma(p+1)*gamma(q+1) / gamma(p+q+4)


###########################
###########################

# 5. qncubed3: Complex Analysis: Double Keyhole Contour
# https://www.youtube.com/watch?v=LT3jpvWMH2s
# - Contour: double keyhole (not a dogbone contour);
integrate(function(x) 1 / (x*sqrt(x^2 - 1)), lower=1, upper=Inf)
pi / 2

### Gen 1:
k = 3
integrate(function(x) 1 / (x * (x^2 - 1)^(1/k)), lower=1, upper=Inf)
pi/2 / sin(pi/k)

###
k = sqrt(5)
integrate(function(x) 1 / (x * (x^2 - 1)^(1/k)), lower=1, upper=Inf)
pi/2 / sin(pi/k)

# Res = 2i*pi * exp(-1i*pi/k)


### Gen 2:
n = 3
k = 2
integrate(function(x) 1 / (x * (x^n - 1)^(1/k)), lower=1, upper=Inf)
pi/n / sin(pi/k)


###
n = sqrt(5)
k = sqrt(2)
integrate(function(x) 1 / (x * (x^n - 1)^(1/k)), lower=1, upper=Inf)
pi/n / sin(pi/k)

