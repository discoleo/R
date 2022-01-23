########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Modular Arithmetic
### Derivation & Experiments
###
### draft v.0.1d




######################

### Helper Functions

source("Polynomials.Helper.R")
# - is automatically loaded in: Polynomials.Helper.R;
#   source("Polynomials.Helper.Factorize.R")


#######################
#######################


### Inv

###
inv.mod(5, 13)
r = 13 %% 5
13 - (13 * inv.mod(r, 5) - 1) / 5

###
inv.mod(5, 17)
r = 17 %% 5
17 - (17 * inv.mod(r, 5) - 1) / 5

###
inv.mod(5, 19)
r = 19 %% 5
19 - (19 * inv.mod(r, 5) - 1) / 5

###
x = 7
p = 19; # gcd(x, p) == 1;
inv.mod(x, p)
r = p %% x
p - (p * inv.mod(r, x) - 1) / x



### SQRT

### 2^2 => sqrt(10) = 1/2 (mod 13)

# 1 (mod 13)
(4 * -3) %% 13
(10 - inv.mod(4, 13))
# x = sqrt(10) (mod 13)
x = inv.mod(2, 13)
x = c(x, 13 - x)
(x*x - 10) %% 13
print(x);

# * 3 => sqrt(12)
(10 * 9) %% 13 # 12 (mod 13)
x2 = (3 * x) %% 13
(x2*x2 - 12) %% 13
print(x2);

# * 5 => sqrt(3)
(10 * 25) %% 13 # 3 (mod 13)
x2 = (5 * x) %% 13
(x2*x2 - 3) %% 13
print(x2);

# * 7 => sqrt(9) # trivial
(10 * 49) %% 13 # 3 (mod 13)
x2 = (7 * x) %% 13
(x2*x2 - 9) %% 13
print(x2);


##########

### 2^2 => sqrt(8) = 1/2 (mod 31)
# 2^4 also possible (mod 31);

p = 31
# 1 (mod p)
(4 * 8) %% p
(8 - inv.mod(4, p))
# x = sqrt(8) (mod p)
x = inv.mod(2, p)
x = c(x, p - x)
(x*x - 8) %% p
print(x);


###
p = 67
#
r =  p %% 4;
r = if(r == 1) 1 else -1;
x2 = (p - r) / 4;
if(r > 0) x2 = p - x2;
(x2 - inv.mod(4, p))
# x = sqrt(x2) (mod p)
x = inv.mod(2, p)
x = c(x, p - x)
(x*x - x2) %% p
print(x2); print(x);

