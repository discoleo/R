
###
### Multi-Clustering
###
### Leonard Mada
###
### draft v.0.1a


### References

# based on:
# Co-manifold learning with missing data - Eric Chi
# https://www.youtube.com/watch?v=TUjgoRswNTc


###############

### Theory

### Cost-Function:
# F(u) = a * sum(D(x[i], u[i])) + g * sum(D(u[i], u[j]));
# Task: min(F)
# where:
# x[i] = Data point i;
# u[i] = Cluster i;
# D(x[i], u[i]) = distance (metric) between point x[i] and "cluster" u[i];
# D(u[i], u[j]) = distance (metric) between cluster u[i] and cluster u[j];
# a, g = coefficients;

### Points in Cluster
### D(x[i], u[i])
# - could be a polynomial combination of Lp norms;
# D(x[i], u[i]) = b2*(x[i] - u[i])^2 + b1*abs(x[i] - u[i]);

### Clusters
### D(u[i], u[j])
# - Regularization terms to merge clusters;
# D(u[i], u[j]) = sum( g[dim, 2] * w[ij] * L2(u[i], u[j]) + g[dim, 1] * w[ij] * L1(u[i], u[j])) )
# - could be a polynomial combination of Lp norms;
# dim = 2D, 3D, 4D;
# w[ij] = weighting / normalisation;


#############
### Types ###
#############

### Bi-Clustering: 2D
F(u) = a * sum(D(x[i], u[j])) + g[1] * sum(Dx(u[i], u[j])) + g[2] * sum(Dy(u[i], u[j]));


### 3D Clustering
F(u) = a * sum(D(x[i], u[j])) + g[1] * sum(Dx(u[i], u[j])) + g[2] * sum(Dy(u[i], u[j])) + g[3] * sum(Dz(u[i], u[j]));



### 4D Clustering
F(u) = a * sum(D(x[i], u[j])) + # individual points
	g[1] * sum(Dx(u[i], u[j])) + g[2] * sum(Dy(u[i], u[j])) + g[3] * sum(Dz(u[i], u[j])) + # (x,y,z)-components
	g[4] * sum(Dt(u[i], u[j])); # time-component / dynamical component;



