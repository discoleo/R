
### Melting Temperature
# https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/documents/367/000/meltingtemp1.pdf


# Note: dS = actual_dS / 1000;
Tm.data = data.frame(
	nn = c("AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG", "CA", "CT", "CC", "CG",
		"GA", "GT", "GC", "GG"),
	dH = c(-9.10, -8.60, -6.50, -7.80, -6.00, -9.10, -5.60, -5.80, -5.80, -7.80,
		-11.00, -11.90, -5.60, -6.50, -11.10, -11.00),
	dS = c(-0.0240, -0.0239, -0.0173, -0.0208, -0.0169, -0.0240, -0.0135, -0.0129,
		-0.0129, -0.0208, -0.0266, -0.0278, -0.0135, -0.0173, -0.0267, -0.0266) );

split.s2 = function(x) {
	splitNc = function(x) {
		tmp = strsplit(x, "(?<=..)", perl = TRUE)[[1]];
		LEN = length(tmp);
		if(nchar(tmp[LEN]) == 1) tmp = tmp[ - LEN];
		return(tmp);
	}
	s1 = splitNc(x);
	s2 = substr(x, 2, nchar(x));
	s2 = splitNc(s2);
	sA = c(s1, s2);
	return(sA);
}

# nc = Oligonucleotide concentration;
Tm = function(x, nc = 0.5, cNa = 0.05) {
	s2 = split.s2(x);
	id = match(s2, Tm.data$nn);
	dH = sum(Tm.data$dH[id]);
	dS = sum(Tm.data$dS[id]);
	# Scaled by 1/1000;
	A = -0.0108;  # dS for Initiation
	R = 0.001987; # Gas Constant
	# 16.6 / log(10) = 7.2093
	cNc = nc * 1E-6; # in mol / l;
	div = A + dS + R*log(cNc/4);
	Tm  = dH / div - 273.15 + 7.2093 * log(cNa);
	return(Tm);
}
# Seq of Nucleotides
Tm.nnSeq = function(x, nc = 0.5, cNa = 0.05) {
	LEN = length(x);
	if(LEN == 0) return(numeric(0));
	if(LEN == 1) return(NA);
	s2 = paste0(x[-LEN], x[-1]);
	id = match(s2, Tm.data$nn);
	dH = sum(Tm.data$dH[id]);
	dS = sum(Tm.data$dS[id]);
	# Scaled by 1/1000;
	A = -0.0108;  # dS for Initiation
	R = 0.001987; # Gas Constant
	# 16.6 / log(10) = 7.2093
	cNc = nc * 1E-6; # in mol / l;
	div = A + dS + R*log(cNc/4);
	Tm  = dH / div - 273.15 + 7.2093 * log(cNa);
	return(Tm);
}

### Simulations
simTm = function(n, iter = 2000, probs = c(1/4,1/4,1/4,1/4),
		nn.seq = c("A","T","G","C")) {
	replicate(iter, {
		nn = sample(nn.seq, n, replace = TRUE, prob = probs);
		Tm.nnSeq(nn);
	});
}
### Sequence = exactly same composition;
simTm.nnSeq = function(nn.seq, iter = 2000) {
	replicate(iter, {
		LEN = length(nn.seq);
		nn  = sample(nn.seq, LEN);
		Tm.nnSeq(nn);
	});
}


### Tests
Tm("AAAAACCCCCGGGGGTTTTT")
Tm.nnSeq(strsplit("AAAAACCCCCGGGGGTTTTT", "", fixed=TRUE)[[1]])
# 69.67

### Simulations
nn = rep(c("A","T","G","C"), c(4,4,4,4))
tm = simTm.nnSeq(nn)
hist(tm, breaks = 20)
summary(tm)

