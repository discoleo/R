

# Polynomial Systems

***Leonard Mada***


## A.) Symmetric Systems

* if (x,y,z) is a solution, so is every permutation;

***Examples***
~~~
x^n + y^n + z^n = R1
x*y + x*z + y*z = R2
x*y*z = R3
~~~


## B.) Hetero-Symmetric Systems

* if (x,y,z) is a solution, so are all ***cyclic*** permutations of it;

### B.1.) Simple Ht-Symmetric

***Examples***
~~~
x1^n*x2^m + x2^n*x3^m + x3^n*x4^m + x4^n*x1^m = R1
x1*x2 + x2*x3 + x3*x4 + x4*x1 = R2
x1*x2*x3 + x2*x3*x4 + x3*x4*x1 + x4*x1*x2 = R3
x1*x2*x3*x4 = R4
~~~

- (usually) does **NOT** permit trivial solutions, like: x1 = x2 = x3 = x4;


### B.2.) Special Ht-Symmetric

***Examples***
~~~
x1^n + b*x2 = R
x2^n + b*x3 = R
x3^n + b*x4 = R
x4^n + b*x1 = R
~~~

- includes also the trivial solution: x1 = x2 = x3 = x4;
