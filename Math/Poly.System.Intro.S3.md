
# Systems of Polynomial Equations with 3 variables and Cyclic Symmetry

***Leonard Mada***,
***...***

We will present a special technique to solve systems of polynomial equations with 3 variables and cyclic symmetry.

A polynomial function F(x,y,z) has cyclic symmetry if the following condition is satisfied:
```
F(x,y,z) = F(y,z,x) = F(z,x,y);
```

We will extend this concept to systems of polynomial equations.

## Classification:
There are 2 major types of systems of polynomial equations with cyclic symmetry: with 3 independent equations or with 1 single type of equation.

**Case 1: System with 3 independent equations**
```
F1(x, y, z) = 0
F2(x, y, z) = 0
F3(x, y, z) = 0
```
Where F1, F2, F3 are 3 distinct polynomial functions with cyclic symmetry.

**Case 2: System with a single type of equation**
```
F(x, y, z) = 0
F(y, z, x) = 0
F(z, x, y) = 0
```
Where the polynomial function F(x,y,z) does not have cyclic symmetry. The system is obtained by cycling the equation 3 times.

Notice that if (x, y, z) is a solution to any of these systems, then the cyclic permutations are also solutions to the original system.

Let S, E2, E3 be the elementary symmetric polynomials of 3 variables:
```
S = x + y + z;
E2 = x*y + y*z + z*x;
E3 = x*y*z;
```

Therefore, S, E2 and E3 are invariant under any cyclic permutation of the original solutions. We will show how to transform the original system into a system in the variables S, E2, E3. After solving this system, one retrieves the solution (x,y,z) by solving the simple polynomial of order 3:
```
x + y + z = S;
x*y + y*z + z*x = E2;
x*y*z = E3;
```

Regarding systems of type 1, there are 2 sub-variants of such systems.

**Case 1.a:** If all of F1, F2, and F3 are symmetric polynomials, then the corresponding system can be easily decomposed into 2 simpler systems:
```
P1(S, E2, E3) = 0
P2(S, E2, E3) = 0
P3(S, E2, E3) = 0,
```
where S, E2, E3 are the elementary symmetric polynomials and P1, P2 and P3 are the corresponding decomposition of F1, F2 and F3. After solving this system, one retrieves the roots by solving the simple polynomial of order 3.

**Case 1.b:** If at least 1 of F1, F2 or F3 is not fully symmetric (but still with cyclic symmetry), then the decomposition of the original system is more involved. The details will be provided further below.

TODO

Ref. GitHub:
- see also: https://github.com/discoleo/R/tree/master/Math/Poly.System.md;

**Type 3 Equations:**
- Fully Symmetric: Poly.System.Symmetric.P3_Full.R
- Hetero-Symmetric Poly.System.Hetero.Symmetric.S3.Mixed.R

**Type 1 Equation:**
- Poly.System.Hetero.Symmetric.S3.R

**Other types:**
- Poly.System.Hetero.Symmetric.S3.Composite.R

