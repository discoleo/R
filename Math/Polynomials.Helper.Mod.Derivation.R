########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Modular Arithmetic
### Derivation & Experiments
###
### draft v.0.1a




######################

### Helper Functions

source("Polynomials.Helper.R")
# - is automatically loaded in: Polynomials.Helper.R;
#   source("Polynomials.Helper.Factorize.R")


#######################
#######################


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

