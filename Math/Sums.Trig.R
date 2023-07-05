

### sum( cos(2*n*x) / n )
# Maths 505:  My take on this on wonderful infinite series from @drpeyam
# https://www.youtube.com/watch?v=mqPTvELJPM0

n = 20000
id = seq(n)

###
x = sqrt(2)
sum(cos(2*x*id)/id)
- log(2*abs(sin(x)))

###
x = sqrt(3)
sum((-1)^id * cos(2*x*id)/id)
- log(2*abs(cos(x)))

### Ex: H(3*n) - H(n)
# 1 + 1/2 - 2/3 + ...
x = pi/3
sum(cos(2*x*id)/id)
- log(2*abs(sin(x)))
#
id = seq(30000)
sum(c(1,1,-2)/id)
log(3)


########################
########################

### sum(cos(n) / n)
# Michael Penn: An infinite cosine sum.
# https://www.youtube.com/watch?v=R_Uf78si8jk

id = seq(60000)
sum(cos(id) / id)
- log(2 - 2*cos(1)) / 2
- log(2*sin(1/2))


