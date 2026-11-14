

n = 30;
r = 10;
x0 = 5; y0 = 2;
#
phi = (2*pi)*runif(n);
dri =  runif(n, -1,1) / 10 * r;
x = x0 + (r+dri) * sin(phi);
y = y0 + (r+dri) * cos(phi);

# Plot:
plot(x, y)

### Circular Regression

#
xm = mean(x); ym = mean(y);
dx = x - xm; dy = y - ym;

# Eq 1:
sum((x-x0)^2 + (y-y0)^2 - r^2) # = 0
r^2 - sum((x-x0)^2 + (y-y0)^2) / n # = 0

# Eq 2 & 3:
sum(((x-x0)^2 + (y-y0)^2 - r^2) * (x-x0)) # = 0
sum(((x-x0)^2 + (y-y0)^2 - r^2) * (y-y0)) # = 0

# Substituting Eq 1 in Eqs 2&3:
sum(((x-x0)^2 + (y-y0)^2 - r^2) * x) # = 0
sum(((x-x0)^2 + (y-y0)^2 - r^2) * y) # = 0
# =>
sum((x-x0)^2*x + (y-y0)^2*x - xm*(x-x0)^2 - xm*(y-y0)^2) # = 0
sum((x-x0)^2*y + (y-y0)^2*y - ym*(x-x0)^2 - ym*(y-y0)^2) # = 0

# =>
sum(x^2*(x-xm) - 2*x*(x-xm)*x0 + (x-xm)*x0^2 + (x-xm)*y^2 - 2*(x-xm)*y*y0 + (x-xm)*y0^2) # = 0
sum(x^2*(y-ym) - 2*x*(y-ym)*x0 + (y-ym)*x0^2 + (y-ym)*y^2 - 2*(y-ym)*y*y0 + (y-ym)*y0^2) # = 0

# =>
sum(x^2*dx + dx*y^2 - 2*x*dx*x0 - 2*dx*y*y0) # = 0
sum(x^2*dy + dy*y^2 - 2*x*dy*x0 - 2*dy*y*y0) # = 0
# =>
sum(x^2*dx + dx*y^2 - 2*x*dx*x0) * sum(dy*y) - sum(x^2*dy + dy*y^2 - 2*x*dy*x0) * sum(dx*y) # = 0
#
div = 2*(sum(x*dx)* sum(dy*y) - sum(x*dy) * sum(dx*y));
(sum(x^2*dx + dx*y^2) * sum(dy*y) - sum(x^2*dy + dy*y^2) * sum(dx*y)) / div;
x0;
#
div = 2*(sum(dx*y)* sum(dy*x) - sum(dy*y) * sum(dx*x));
(sum(x^2*dx + dx*y^2) * sum(dy*x) - sum(x^2*dy + dy*y^2) * sum(dx*x)) / div;
y0;

# TODO: Solve remaining system;

