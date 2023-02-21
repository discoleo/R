

### Integrals: Log(Trig)


Catalan = 0.915965594177219015054603514;
# Note:
# Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)


#################

### Michael Penn: Can you guess the trick for this integral?
# https://www.youtube.com/watch?v=8R0MiRYmjbk
# Intermediary:
#
integrate(function(x) log(tan(x)), 0, pi/4)
# ==
- Catalan
