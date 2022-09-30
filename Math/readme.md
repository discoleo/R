# PolynomOS: a Revolution in Mathematics

## Math Tools & Problems

1.) Polynomials
- Multivariable polynomials in R;
- Generalization of Cardano's formula & exact solutions for higher order polynomials;
- Polynomial systems;

2.) Exact Integrals
- Polynomial Fractions;

3.) Double Integrals
- e.g. with result = Gamma(1/n)^2
- relations of the Gamma function: G(1/2n) = f(G(1/n))

4.) ODEs
- Linear ODEs with Polynomial coefficients;
- Non-linear ODEs;

## Polynomials

### A. Class 1 Polynomials
- P[n]: roots with radicals of order n;

A.1 Generalization of Cardano's formula:
- see file: Polynomials.CardanGeneralisation.R;
- proper way of Cardano's formula: based on roots of unity & independent of coefficients;
- these polynomials are a very special case of Class 1 polynomials;


A.2 Complete Class 1
- see files: Polynomials.Class1.Formulas.R, Polynomials.Class1.R;

### B. Class 2 Polynomials
- Polynomial compositions;

B.1. Polynomial of Order 4
- proper & complete way to solve polynomials of Order 4 (with all 4 roots), see file: Polynomials.P4.R;
- based on the C2-Decomposition, see: Poly.System.S4.C2.R and related files;

B.2 Based on Roots of Unity
- Special case of Class 2: P[n] based on roots of unity of order (n+1), see: Polynomials.Class2.Formulas.R;

B.3. Special Compositions
- based on Hetero-Symmetric Polynomial Systems (cyclic permutations);
- see e.g. Poly.System.Hetero.Symmetric.R, Poly.System.Hetero.Symmetric.S3.R, Poly.System.Hetero.Symmetric.S4.R;

### C. Class 3 Polynomials
- Roots based on cos(2*pi/n);
