

source("PCR.R")


### Tests
Tm("AAAAACCCCCGGGGGTTTTT")
Tm.nnSeq(strsplit("AAAAACCCCCGGGGGTTTTT", "", fixed=TRUE)[[1]])
# 69.67


### Example 2:
Tm("AAAACCCCGGGGTTTT")


### Various Tm: 4 G, 4 C, 4 A, 4 T
nn = rep(c("A","T","G","C"), c(4,4,4,4))
r = find.simTm(nn, 55, iter = 10000, tol = 0.02)
print(r); print(Tm.nnSeq(r));
r = find.simTm(nn, 40, iter = 10000, tol = 0.1)
print(r); print(Tm.nnSeq(r));

Tm("AAAACCCCGGGGTTTT") # 61.5 C
Tm("TGGCATGCTTCCAGAA") # 60 C
Tm("ATGGGGCCTTCAATCA") # 60 C
Tm("CCAGGTTCGGCATATA") # 55 C
Tm("AGTCGACATTCAGCTG") # 52 C
Tm("ACCTCTGTAGACTAGG") # 42 C
Tm("CCAGGAGTCTACTAGT") # 42 C
Tm("GACTCTGAGTACACTG") # 40 C
Tm("GTAGTCCGTACAGTAC") # 40 C

### Ionic Concentration

nn = "CCAGGTTCGGCATATA";
Tm(nn, cNa = 0.100) # 60.0 C
Tm(nn, cNa = 0.075) # 57.9 C
Tm(nn, cNa = 0.050) # 55.0 C # Standard conc.
Tm(nn, cNa = 0.025) # 50.0 C


###################
### Simulations ###

###
nn = rep(c("A","T","G","C"), c(4,4,4,4))
tm = simTm.nnSeq(nn)
hist(tm, breaks = 20)
summary(tm)


###
# png("PCR.Tm.16nn.png")
nn = rep(c("A","T","G","C"), c(4,4,4,4))
tm = simTm.nnSeq(nn, iter = 24000)
hist(tm, breaks = 20, main = "Histogram Tm", xlab="Tm (\uBA C)", cex=2)
summary(tm)
# dev.off()


# png("PCR.Tm.16nn.png")
par.old = par(cex.main = 2, cex.lab = 2)
# hist(tm, breaks = 20, main = "Histogram Tm", xlab="Tm (\uBA C)", cex=2)
hist(tm, breaks = 20, main = "Histogram Tm", xlab=NULL, cex=2)
par(par.old)
# dev.off()

#######################

### Primers:

n = 200
x = sample(c("G","C","A","T"), n, TRUE)
xinv = complement.nn(rev(x))
x = paste(x, collapse = "")

p = find.primer(x)
print(p$Match)


############################

### Examples from Literature

### Murine EGFR:
# - P1: has probably an additional sequence;
Tm("ACCGCTAGCGCCGCCACCATGCGACCCTCAGGGACCGC")
Tm("ACCAAGCTTTCATGCTCCAATAAACTCACT")

# Set 2:
Tm("AGTTGTTATCAGTAAGGGAGCTGCA")
Tm("ACCGAAAATCTGTGGGAAGTCTTGT")

# PCR: 95ºC 3 min, [95ºC 30 s, 58ºC 30 s, 72ºC 40 s] x 38 cycles
# 72ºC 3 min, 4ºC 5 min;


### TP53

Tm("GGTTAAACCCAGCTTGACCA")
Tm("GGAGGCAGAGACAGTTGGAG")

# PCR: 94ºC 1 min, [94ºC 30sec, 60ºC 30 s, 72ºC 1 min] x 35 cycles
# 72ºC 10 min, 4ºC 5 min;


######################
######################


### PCR / qPCR

# png(file = "PCR.Quantity.Cycles.png")

# Exponential growth, then levelling off;
plot.quantity()

# dev.off()
