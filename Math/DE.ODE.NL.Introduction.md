
# Non-Linear (NL) ODEs

**Leonard Mada**

## Introduction

This material presents various types of Non-Linear (NL) ODEs. The focus is on NL ODEs with polynomial coefficients. Both homogenous and inhomogeneous variants are covered and an explanation will be provided how to generate both types.

### Motto:
: The new AI proofs in Maths look only more and more ridiculous!

## Types

The NL ODEs will be organized by the leading Order term.

### Type: y * d2y

There are various ways to generate this type of NL ODEs. There are 2 broad ways to generate such ODEs:
- Simple ways: using the power of a function;
- Coupled ways: Coupling the power of a function with an additional function;

### Simple variants

1. Simple Power of Log
```
y = B(x) * Log(P(x))^n
```

2. Power of Log w. Exp
```
y = B(x) * Log(P1(x))^n * Exp(P2(x))
```

For examples, see file: DE.ODE.NL.Log.PowExp.R;


3. Power of Log w. Radical
```
y = B(x) * Log(P1(x))^n1 * P2(x)^(1/n2)
```

For examples, see file: DE.ODE.NL.Log.PowRadical.R;


4. Power of Atan
```
y = B(x) * Atan(P(x))^n
```

For examples, see file: DE.ODE.NL.Atan.Pow.R;
TODO: generalize power;


5. Power of Atan w. Exp
```
y = B(x) * Atan(P1(x))^n * Exp(P2(x))
```


6. Power of Atan w. Radical
```
y = B(x) * Atan(P1(x))^n * P2(x)^r
```


7. Exp of Radical
```
y = C(x) * Exp(B1(x) * P(x)^r + B0(x))
```

For examples, see file: DE.ODE.NL.Exp.Radicals.R;


8. Exp of Log^n
```
y = C(x) * Exp(Log(P(x))^r)
```

For examples, see file: DE.ODE.NL.Exo.LogPow.R;


#### Coupled Variants

6. Power of Log: Coupled w. Hyperbolic Functions
```
# Note: the same P(x);
y = B(x) * Log(Cosh(P(x)))^n * Cosh(P(x))
```

7. Power of Atan: Coupled w. Hyperbolic Functions
```
# Note: the same P(x);
y = B(x) * Atan(Exp(P(x)))^2 * Cosh(P(x))
```

For examples, see file: DE.ODE.NL.Atan.CoupledExp.R; not yet generalized for higher powers.


8. Power of Log: Coupled w. Sqrt
```
# Note: the same P(x);
y = B(x) * Log(Sqrt(P(x)^2 + b) - P(x))^n * Sqrt(P(x)^2 + b)
```

For examples, see file: DE.ODE.NL.LogSqrt.Pow.R;


**TODO:** more variants;
