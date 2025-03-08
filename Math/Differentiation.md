
# Symbolic Differentiation in R

It is possible to differentiate an expression in R using the function **D**.

## Symbolic Differentiation:

**Examples:**

```
D(expression(k1 * Y * (13 - Y - c1*N) / 13), "Y")

# Complicated expression:
D(expression((x^2 - 3*x) * atan(2*exp(-3*x)+1)), "x")
```

**Advanced Expression:**
```
D.expr = function(x, by) {
	y = substitute(expression(x));
	y = y[[2]];
	dx = lapply(by, function(nm) {
		D(y, nm);
	})
	return(dx);
}

### Test
D.expr(x^3 - 4*x*y + 2*y^2 + 5*sin(cos(3*x)), c("x", "y"))
```

**Note:**
- Result is not simplified;


## Substitution of Parameter

### 1. On Result of D
```
x = D(expression(k1 * Y * (13 - Y - c1*N) / 13), "Y")
params = list(k1 = 2, c1 = 3)
do.call(substitute, list(x, params))
```

### 2. On Initial Expression

```
params = list(k1 = 2, c1 = 3)
x = substitute(D(expression(k1 * Y * (13 - Y - c1*N) / 13), "Y"), params)
eval(x)
```

### 3. Evaluate Expression
Requires all parameters & variables to be assigned a value.


