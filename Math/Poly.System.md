

# Polynomial Systems

***Leonard Mada***

Various classes of polynomial systems can be decomposed into systems of lower order.

These systems of lower order provide a means to solve higher order polynomials using lower order polynomials which are entangled into a very complicated way. The details of the entanglement emerge during the procedure to solve the lower order polynomial system.

The 2 major classes of such systems are:
1. Fully symmetric systems;
2. Hetero-symmetric systems: these are invariant under a cyclic permutation;


## A.) Symmetric Systems

* If (x,y,z) is a solution, so is every permutation of it;
* see files:
  * Poly.System.Symmetric.P3_Full.R;
  * Poly.System.Symmetric.S3.Multiple.R;

***Examples***
~~~
# 3 Variables
x^n + y^n + z^n = R1
x*y + x*z + y*z = R2
x*y*z = R3
~~~

In this particular example, the system can be decomposed into 2 sequential polynomials: P\[n] o P\[3], where P\[n] is a polynomial of order n, and P\[3] is a polynomial of order 3 whose coefficients are functions of the roots of P\[n].

Such systems are a very special case of the systems described in section \[B.1] (and are much easier to solve).

The composition operator "o" is an inverse composition:
```
P[n] o P[3] = P[n](P[3]);
```

**Note:** This variant can be generalised to any number of variables and any type of fully symmetric equations. The particular example above is just very easy to solve.

---

## B.) Hetero-Symmetric Systems

* If (x,y,z) is a solution, so are all ***cyclic*** permutations of it;

There are 2 main variants of this type.

### B.1.) Simple Ht-Symmetric

***Examples***
~~~
x1^n*x2^m + x2^n*x3^m + x3^n*x4^m + x4^n*x1^m = R1
x1*x2 + x2*x3 + x3*x4 + x4*x1 = R2
x1*x2*x3 + x2*x3*x4 + x3*x4*x1 + x4*x1*x2 = R3
x1*x2*x3*x4 = R4
~~~

- (usually) does **NOT** permit trivial solutions, like: x1 = x2 = x3 = x4;
- Some equations may be fully symmetric: but this is not a requirement;

The generalized example:
- Let P be the class of polynomial functions which are invariant under the cyclic permutation of their variables;
- Let P1, ..., Pn be n such functions, each with n variables: x1, ..., xn;
- Let Sys(P, n) be the system formed by these n equations: Pi(x1, ..., xn) = 0;
- This system can be decomposed into 2 subsystems of lower order;

For more examples, see also:
- Poly.System.Hetero.Symmetric.S3.Mixed.R;


### B.2.) Special Ht-Symmetric

***Examples***
~~~
x1^n + b*x2 = R
x2^n + b*x3 = R
x3^n + b*x4 = R
x4^n + b*x1 = R
~~~

This system has 2 trivial solutions:
- Sys\[2, n]: x1 = x3 & x2 = x4, becoming a system of 2 variables;
- Sys\[1, n]: trivial solution: x1 = x2 = x3 = x4, which is also a trivial solution to Sys\[2,n];

This particular system can be actually decomposed into several independent polynomials/systems:
```
- P[n] * P[n^2-n] * P[n^4 - n^2];
- P[n^2 - n] can be solved by solving 2 lower order systems: P[(n^2-n)/2] o P[2];
- P[n^4 - n^2] can be solved by solving 2 lower order systems: P[(n^4-n^2)/4] o P[4];
```

The generalized example:
- Let F be a polynomial function with n variables;
- Let Sys(F, n) be the system formed from the following equations:
```
F(x1, ..., xn) = 0;
F(x2, ..., xn, x1) = 0;
...
F(xn, ..., x[n-1]) = 0;
```

For more details and examples, see also:
- Poly.System.Hetero.Symmetric.S3.Theory.R: but only minimalistic information;
- Poly.System.Hetero.Symmetric.S3.R: the basic examples for the S[3] system;
- most files may need extensive cleanup!


## Solutions

General strategy to solve these systems:

- Observation 1: If (x1, ..., xn) is a solution to such a system, then every cyclic permutation is also a solution.

- Observation 2: The elementary polynomials formed by (x1, ..., xn) are invariant under a cyclic permutation;
```
# Elementary-System
S = x1 + ... + xn;
E2 = x1*x2 + x1*x3 + ... + x[n-1]*xn;
...
En = x1*x2*...*xn;

# Note:
# - S Notation: introduced to distinguish it easily from
#   the other elementary polynomials;
# - Extensively used in the solutions: easier to write;
```

All these polynomials are invariant under the cyclic permutation of the solution.

The main strategy is to transform the original system into a system in {S, E2, ..., En} and solve this lower order system.

Step 2: solve the Elementary system;
- this is usually a polynomial of Order n = number of variables;

Note: Systems of Type \[B.2] are in general easier to transform into the lower order system than systems of type \[B.1].

There are many examples on GitHub:
- see most files named Poly.System.S[n].Ht...R, where n is the number of variables;
