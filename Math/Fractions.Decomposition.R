


### Fraction Decopositions
### Other Fractions

# - for 1 / (x^n + 1):
#   see file Integrals.Fractions.Unity.R;
# - for 1 / CardanPoly[n]:
#   see file Integrals.Fractions.CardanPoly.R;


###################

### Basic Fractions

### 1 / (x^n * (x + k))

### n = 2
k = sqrt(3)
x = pi/exp(1)
#
1/(x^2 * (x+k))
1/k * 1/x^2 - 1/k^2 * 1/x + 1/k^2 * 1/(x+k)

# Derivation:
1/(x^2 * (x+k))
1/k * 1/x * (1/x - 1/(x+k))
1/k * (1/x^2 - 1/(x*(x+k)))
1/k * 1/x^2 - 1/k^2 * (1/x - 1/(x+k))
1/k * 1/x^2 - 1/k^2 * 1/x + 1/k^2 * 1/(x+k)


### n = 3
k = sqrt(3)
x = pi/exp(1)
#
1/(x^3 * (x+k))
1/k * 1/x^3 - 1/k^2 * 1/x^2 + 1/k^3 * 1/x - 1/k^3 * 1/(x+k)

#
n = 4
1/(x^n * (x+k))
1/k * 1/x * (1/x^n - (-1)^n/k^n) / (1/x + 1/k) + (-1)^n/k^n * 1/(x+k)

