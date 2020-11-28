
setwd("C:\\Users\\Leo Mada\\Desktop\\DB\\Geo")

### Various Integrals
### Unity & Variants:
### 1/(x^n + 1)
### but NOT correctly solved!


i <- 1:50000
k <- 3

x.f <- function(x) prod((1 + 1/i^k)/(1+x^i))

### sequential
x.f.f <- function(x) {
  rez <- double(length(x))
  for(i in 1:length(x)) {
    rez[i] <- x.f(x[i])
  }
  return(rez)
}
### parallel
x.f.f <- function(x) {
  # rez <- double(length(x))
  rez = sapply(1:length(x), function(i) {
    # rez[i] <<- x.f(x[i])
	x.f(x[i])
  })
  return(rez)
}

x <- integrate(x.f.f, lower=0, upper=1)


x$value * 2


exp(x$value) * 2

F(n) = Integral[x = 0, 1] Prod[i = 1, n] (1/(1+x^i)) dx

F(1) = ln(2)
F(2) = 1/2 * [ln(2)/2 + arctan(1)] = ln(2)/4 + pi/8

i <- 1:2

x.f <- function(x) prod(1/(1+x^i))

x <- integrate(x.f.f, lower=0, upper=1)
x <- x$value
log(2)/4 + pi/8 - x

############################



int.f <- function(n) {
id = 1

for(n.pow in n) {
   x.f <- function(x) 1/(1 + x^n.pow)
   rez.df[id, "rez"] <<- integrate(x.f, lower=0, upper=1)$value
   id <- id + 1
}
}

n <- 2^(1:10)
rez.df <- data.frame(n=n, rez=NA)

int.f(n)
rez.df

### Even
n <- 2*(2*(0:30)+1)
rez.df <- data.frame(n=n, rez=NA)

int.f(n)
rez.df
plot(rez.df)

### Odd
n <- 2*(1:30)+1
rez.df <- data.frame(n=n, rez=NA)

int.f(n)
rez.df
plot(rez.df)

############################


root.f <- function(pow, beta, k) {
print(paste("Polynomial: x^", pow, " - ", beta, " * x - ", k, " = 0", sep=""))

r.n <- k^(1/pow)
### Approximate solution
print("Approximate solution:")
r.approx <- beta * r.n / (r.n^(pow - 1) * pow - beta) + r.n
r <- r.approx
err <- r^pow - beta * r - k

print(paste("Root = ", r, ", Error = ", err, sep=""))

### Accurate solution
r <- k^(1/pow)
err <- r^pow - r - k
it <- 1
while(TRUE) {
  it <- it + 1
  r <- (r * beta + k)^(1/pow)
  err <- r^pow - beta * r - k

  print(paste("It = ", it, ", Root = ", r, ", Error = ", err, sep=""))
  if(abs(err) <= 1E-10) {
    break
  }
}

return(r)
}

root.f(5, 1, 1)
root.f(5.5, 1, 7)
root.f(5.5, -1, 7)

root.f(5, 2, 5)
root.f(5, -2, 5)

###############################

g.f <- function(x, nn) {
	r <- gamma(x)
	r <- 1/r^nn
	return(r)
}

g.f.i <- function(pow) {
	integrate(g.f, 0.0001, Inf, pow, subdivisions=10000)
}

g.f.i(1)

g.div.f <- function(n, div=1.5) {
	r <- integrate(g.f, div, Inf, n, subdivisions=10000)$value/
		integrate(g.f, 0.00001, div, n, subdivisions=10000)$value
		return(r)
}

g.div.f(1)
g.div.f(1, exp(1))
[1] 0.2657412

###

# infinite Product
i <- 1:30000

x.f <- function(x, n) prod(1/(1+x^(n*i)))

x.f.f <- function(x, n) {
  rez <- double(length(x))
  for(i in 1:length(x)) {
    rez[i] <- x.f(x[i], n)
  }
  return(rez)
}


x <- integrate(x.f.f, lower=0, upper=Inf, 1, subdivisions=1000)
x


### Solution
xp.f <- function(x, p) prod(1 - x^p)
x.p.f <- function(x, p) {
  rez <- double(length(x))
  for(i in 1:length(x)) {
    rez[i] <- xp.f(x[i], p)
  }
  return(rez)
}
#
prim <- c(1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 91, 97, 101)
#
x.if <- function(upper, primes=prim) {
	x <- integrate(x.p.f, lower=0, upper=upper, primes, subdivisions=4000)$value
	return(x)
}
x <- x.if(1)
x

1/1 - 1/2 - 1/4 + 1/5 - 1/6 + 1/7 + 1/9 - 1/10 - 1/8 + 1/9 + 1/11 - 1/12 + 1/13 - 1/14 - 1/16 + 1/17
[1] 0.4328065

###

x0 <- integrate(x.p.f, lower=0, upper=1/2, prim, subdivisions=1000)$value

plot(tapply(1:50, 1:50, function(id) x/integrate(x.p.f, lower=0, upper=1/id, prim, subdivisions=1000)$value))

###

i <- 1:30000

x.f <- function(x, n) prod((1 + 1/i^n)/(1+x^i))

x.f.f <- function(x, n) {
  rez <- double(length(x))
  for(i in 1:length(x)) {
    rez[i] <- x.f(x[i], n)
  }
  return(rez)
}

x <- integrate(x.f.f, lower=0, upper=1, 2)
i <- 1:30000

x.f <- function(x, n) prod((1 + 1/i^n)/(1+x^i))

x.f.f <- function(x, n) {
  rez <- double(length(x))
  for(i in 1:length(x)) {
    rez[i] <- x.f(x[i], n)
  }
  return(rez)
}

x <- integrate(x.f.f, lower=0, upper=1, 2)
x


### x = All Radicals/Powers
xpow.f <- function(x, p) prod(1 - p*x^p)
x.pow.f <- function(x, p) {
  rez <- double(length(x))
  for(i in 1:length(x)) {
    rez[i] <- xpow.f(x[i], p)
  }
  return(rez)
}
#
powers <- 1:8192
#
x.pow.if <- function(upper, arrPowers=powers) {
	x <- integrate(x.pow.f, lower=0, upper=upper, arrPowers, subdivisions=1024)$value
	return(x)
}
x <- x.pow.if(0.9)
x

plot.pow.f <- function(x, p=powers) {
	rez <- double(length(x))
  for(i in 1:length(x)) {
    rez[i] <- xpow.f(x[i], p)
  }
  return(rez)
}

curve(plot.pow.f, from=0, to =0.8)
