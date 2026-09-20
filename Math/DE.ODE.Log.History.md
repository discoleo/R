
# History of DE.ODE.Log.R
** Leonard Mada**


## History

### Draft v.0.3t:
- [refactor] Moved History to this file;

### Draft v.0.3m - v.0.3s:
- [refactor]
- P(x)^LOG moved to new file: DE.ODE.NL.Log.PowL.R;
- LOG( LOG ) moved to new file: DE.ODE.NL.Log.Log.R;
- Moved ODEs based on Hidden Log to new file: DE.ODE.NL.Log.Other.R;
- Prod( LOG ) moved to new file: DE.ODE.NL.Log.Prod.R;
- Int( LOG ) moved to new file: DE.ODE.NL.Log.Int.R;
- LOG( OTHER(F(x)) ) moved to new file: DE.ODE.NL.Log.Composite.R;

### Draft v.0.3i - v.0.3j:
- Mixed variants: Sum(LOG, EXP)
: y = B1(x)\*log(P1(x)) + B2(x)\*exp(P2(x)) + F0(x);
- Sum(LOG, I( EXP )):
: y = B1(x)\*log(P1(x)) + B2(x)\*exp(P2(x))*I( exp(-P2(x)) ) + F0(x); [v.0.3j]

### Draft v.0.3h:
- Workout of case:
: y = log(P) \* log(log(P));

### Draft v.0.3g - v.0.3g-ex2:
- Automatic generation of simple types of ODEs;
- More examples & checks; [v.0.3g-ex2]

### Draft v.0.3f - v.0.3f-st1:
- Derived from: [TODO full derivation]
: y = log(exp(P1(x)) + P2(x)) + F0(x); [started]
: y = log(log(P(x))) + F0(x);
: y = log(P(x)) \* log(log(P(x))) + F0(x); [v.0.3f-bis]

### Draft v.0.3e:
- Derived from:
: y \* log(x+k) = I(log(x+k) / x) + F0(x);

### Draft v.0.3d:
- Clean-up;

### Draft v.0.3c - v.0.3c-ex2:
- Mixed Log-Exp:
: y = log(P1(x)) \* exp(P2(x)) + F(x);
- [refactor] moved to file: DE.ODE.Mixed.Exp.Log.R;
- More examples; [v.0.3c-ex2]

### Draft v.0.3a - v.0.3b:
- [refactor] moved to new file;
- from:
: y = (x + k1)^(x + k2) + F0(x);
: y * (x + k1)^(x + k2) = F0(x);
- Example:
:(x+k1)*(y - f0)*d2y = (x+k1)*(dy - df0)^2 + (y - f0)^2 + (x+k1)*d2f0*(y - f0);

### Draft v.0.2c:
Clean-up:
- moved section: y = log(P1(x)) \* log(P2(x)) from DE.ODE.Fractions.Lambert.R to "this" file;

### Draft v.0.2b:
- Derived from:
: y \* I(1/log(x + k)) dx = F0(x);

### Draft v.0.2a:
- derived from:
: y = x \* I(1/log(x + k)) dx + F0(x);

### Draft v.0.1d - v.0.1g:
- Derived from: y = x^log(x)
: x^2\*y\*d2y - x^2\*dy^2 + x\*y\*dy - 2\*y^2 = 0;
- Extension: y = (x+k)^log(x+k)
: (x+k)^2\*y\*d2y - (x+k)^2\*dy^2 + (x+k)\*y\*dy - 2\*y^2 = 0;
- Generalization to higher orders:
: y = (x^n+k)^log(x^n+k); [v.0.1f]
- Generalization:
: y = P(x)^log(P(x));
