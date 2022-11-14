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
- proper way of Cardano's formula: based on roots of unity & quasi-independent of coefficients;
- these polynomials are a very special case of Class 1 polynomials;

A.2 Complete Class 1
- see files: Polynomials.Class1.Formulas.R, Polynomials.Class1.R;
- Rationalization of Fractions with radicals, see file: Fractions.Rationalisation.R;


### B. Class 2 Polynomials
- Polynomial compositions, e.g. P[6] = P[2] o P[3] or P[6] = P[3] o P[2];
  - for simple examples of basic compositions, see Polynomials.Derived.P6.R;
- the coefficients of the components can show complex entanglements (see B.3);

B.1. Polynomial of Order 4
- proper & complete way to solve polynomials of Order 4 (with all 4 roots), see file: Polynomials.P4.R;
- based on the C2-Decomposition (P[2] o P[2] o P[3]), see: Poly.System.S4.C2.R and related files;

B.2 Based on Roots of Unity
- Special case of Class 2: P[n] based on roots of unity of order (n+1), see: Polynomials.Class2.Formulas.R;
- Symmetric polynomials are a generalization of the roots of unity sub-class, see: Polynomials.Derived.P6.Symmetric.R;

B.3. Special Compositions
- based on Hetero-Symmetric Polynomial Systems (cyclic permutations);
  > see e.g. Poly.System.Hetero.Symmetric.R, Poly.System.Hetero.Symmetric.S3.R, Poly.System.Hetero.Symmetric.S4.R;


### C. Class 3 Polynomials
- Roots based on cos(2*pi/(2*n+1));


## ODEs

> **A. Linear ODEs**\
> **B. Non-Linear ODEs**\
> **C. Systems of ODEs**

### A. Linear ODEs

Linear ODEs with Polynomial Coefficients are constructed starting from the following functions:
- y = B1 * exp(P1) + B2 * exp(P2) + B0;
- y = B1 * log(P1) + B2 * log(P2) + B0;
- y = B1 * atan(P1) + B2 * atan(P2) + B0;
- Mixed: y = B1 * exp(P1) + B2 * log(P2) + B0 and all other combinations; (TODO)
- y = B1 * sin(P) + B2 * cos(P) + B0; (Note: the same function P(x))
- Mixed Exp-Trig: y = sin(P1) * exp(P2) + B0, see file DE.ODE.Gaussian.R (currently only a simple example);

Note:
- B2(x), B1(x), B0(x), P1(x), P2(x) are polynomials or polynomial fractions of x;
  > Note: the examples are usually based on simple polynomials, as the fractions are much uglier to work with;
