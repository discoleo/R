###################
##
## PCR Tools
##
## Leonard Mada
##
## draft v.0.2a

### PCR Tools
# Functions to:
# - Compute Tm (melting temperature);
# - Estimate Tm in the presence of 1 mismatch;
# - Match optimal primers for PCR;


### Ref:

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

split.s2 = function(x, rm.space = "[\\s]++") {
	splitNc = function(x) {
		tmp = strsplit(x, "(?<=..)", perl = TRUE)[[1]];
		LEN = length(tmp);
		if(nchar(tmp[LEN]) == 1) tmp = tmp[ - LEN];
		return(tmp);
	}
	if(! is.null(rm.space)) {
		x = gsub(rm.space, "", x, perl = TRUE);
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

### Tm with 1 Mismatch
# Note: wild guess, but I do NOT have the time
#  to search the literature;
# npos = location of mismatch;
# type = add penalty to Tm based on type of mismatch:
# - bulky: only if bulky (A/G or G/A);
#   Note: only the base pair in the primer is checked!
# - no: do not add any penalty;
# - always: independent of actual nucleotides;
Tm.mismatch = function(x, npos, type = c("bulky", "no", "always"),
		nc = 0.5, cNa = 0.05, verbose = TRUE) {
	if(npos < 1 || npos > nchar(x)) stop("Invalid position of mismatch!");
	type = match.arg(type);
	x   = strsplit(x, "", TRUE)[[1]];
	LEN = length(x);
	if(npos < 3) {
		Tm = Tm.nnSeq(x[seq(npos + 1, LEN)], nc=nc, cNa=cNa);
		return(Tm);
	}
	if(npos + 2 > LEN) {
		if(verbose) warning("Mismatch is at the 3'-end!");
		Tm = Tm.nnSeq(x[seq(npos - 1)], nc=nc, cNa=cNa);
		return(Tm);
	}
	# dH, dS:
	s2 = paste0(x[-LEN], x[-1]);
	s2 = s2[c(1-npos, -npos)];
	id = match(s2, Tm.data$nn);
	dH = sum(Tm.data$dH[id]);
	dS = sum(Tm.data$dS[id]);
	if((type == "bulky" && x[npos] %in% c("A","G")) ||
		type == "always") {
		nm = c("AG", "GA"); # wild guess;
		id = match(nm, Tm.data$nn);
		mH = sum(Tm.data$dH[id]) / 2;
		mS = sum(Tm.data$dS[id]) / 2;
		dH = dH - mH; dS = dS - mS;
	}
	# Scaled by 1/1000;
	A = -0.0108;  # dS for Initiation
	R = 0.001987; # Gas Constant
	cNc = nc * 1E-6; # in mol / l;
	div = A + dS + R*log(cNc/4);
	Tm  = dH / div - 273.15 + 7.2093 * log(cNa);
	return(Tm);
}

### Complementary Seq:
complement.nn = function(x, rev = FALSE, collapse = NULL) {
	LEN = length(x);
	if(LEN == 0) return(x);
	y = rep("X", LEN);
	y[x == "G"] = "C";
	y[x == "A"] = "T";
	y[x == "C"] = "G";
	y[x == "T"] = "A";
	if(rev) y = rev(y);
	if(! is.null(collapse)) y = paste(y, collapse=collapse);
	return(y);
}


# TODO: NOT yet (fully) implemented;
# x  = String with nucleotide seq;
# Tm = Desired melting temperature;
# is.5p = DNA is in 5' => 3' orientation;
#   => will select primer at 3'-end of DNA;
# skip.nn = can skip first nucleotides;
# type: lead = leading strand; both = both strands;
#  filter = filtered to match both strands;
# TOP = top matches;
# w   = weights used for criteria used to rank matches;
find.primer = function(x, skip.nn = 5, is.5p = TRUE, TOP = 10,
		type = c("filter", "both", "lead", "non-lead"),
		Tm = 55, keep.Tm = 50, w = c(1, 1/8, 1), GCp = 0.5, tol = 0.5) {
	type = match.arg(type);
	# Seq of NN:
	x = strsplit(x, "", fixed = TRUE)[[1]];
	if(! is.5p) x = rev(x);
	if(type == "lead") {
		Rez = find.primer.nnSeq(x, is.5p = TRUE, Tm=Tm, keep.Tm=keep.Tm,
			tol=tol, skip.nn = skip.nn);
		return(Rez);
	}
	xinv = rev(complement.nn(x));
	P2 = find.primer.nnSeq(xinv, is.5p = TRUE, Tm=Tm, keep.Tm=keep.Tm,
			tol=tol, skip.nn = skip.nn);
	if(type == "non-lead") {
		return(P2);
	}
	P1 = find.primer.nnSeq(x, is.5p = TRUE, Tm=Tm, keep.Tm=keep.Tm,
			tol=tol, skip.nn = skip.nn);
	#
	P = list(P1 = P1, P2 = P2);
	if(type == "both") return(P);
	### Filtered:
	# Match Primers:
	P$Match = match.primers(P1, P2, TOP=TOP, Tm=Tm, GCp=GCp, w=w);
	return(P);
}

### Match Primers
# p1, p2 = data.frames with Tm & GCp for the 2 sets of primers;
# TOP = top matches;
# w   = weights applied to criteria used to rank matches;
match.primers = function(p1, p2, TOP = 10, w = c(1, 1/8, 1),
		Tm = 55, GCp = 0.5) {
	nr = nrow(p1);
	if(nr == 0) {
		return(NA);
	}
	flt = sapply(seq(nr), function(id) {
		tmp = p1[id, c("Tm", "GCp")];
		opt = abs(tmp$Tm - p2$Tm) * w[1] +
			abs(tmp$Tm - Tm) * w[2] + abs(p2$Tm - Tm) * w[2] +
			abs(tmp$GCp - GCp) * w[3] + abs(p2$GCp - GCp) * w[3];
	});
	# Order Best Matches:
	id  = order(flt);
	top = min(TOP, length(id));
	id  = id[seq(top)];
	idr = (id - 1) %%  nrow(p2) + 1;
	idc = (id - 1) %/% nrow(p2) + 1;
	Rez = cbind(p1[idc,], P2 = p2[idr,]);
	Rez$dTm   = abs(Rez$Tm - Rez$P2.Tm);
	Rez$Score = flt[id]; rownames(Rez) = NULL;
	return(Rez);
}

### Find Primers
# x = Array of nucleotides;
# Note: for a string of nucleotides, see: find.primer;
find.primer.nnSeq = function(x, is.5p = TRUE, Tm = 55, keep.Tm = 50,
		skip.nn = 5, tol = 0.5, ...) {
	# Leading Strand:
	if(! is.5p) x = rev(x);
	# Non-Leading Strand: TODO
	# Search Complementary DNA-strand;
	# if(is.5p) x = rev(complement.nn(x));
	#
	LEN  = length(x); 
	nLen = 10; nLast = LEN - nLen;
	if(LEN <= 20) {
		stop("Too short! Not yet implemented!");
	}
	nS  = 1;
	# Init Result:
	Rez = list(list(nS=nS, Len = nLen + 1, Tm = list()));
	while(nS <= nLast) {
		id = 1;
		nPosEnd = nS + nLen + id - 1;
		Tmp = Tm.nnSeq(x[nS:nPosEnd], ...);
		Rez[[nS]]$Tm[[id]] = Tmp;
		while(Tmp < Tm) {
			nPosEnd = nPosEnd + 1; id = id + 1;
			if(nPosEnd >= LEN) break;
			Tmp = Tm.nnSeq(x[nS:nPosEnd], ...);
			Rez[[nS]]$Tm[[id]] = Tmp;
		}
		Rez[[nS]]$Tm = unlist(Rez[[nS]]$Tm);
		if(! is.null(keep.Tm)) {
			keep = Rez[[nS]]$Tm >= keep.Tm;
			Rez[[nS]]$Tm  = Rez[[nS]]$Tm[keep];
			Rez[[nS]]$Len = nLen + id + 1 - length(Rez[[nS]]$Tm);
		}
		nS = nS + 1;
		if(nS > skip.nn) break;
		# Init list:
		Rez[[nS]] = list(nS=nS, Len = nLen + 1, Tm = list());
	}
	Rez = as.df.primer(Rez);
	### GC-Content:
	Rez = gc.nn(Rez, x);
	# TODO: actual selection of best candidates;
	return(Rez);
}

### GC-Content:
gc.nn = function(primers, x) {
	Rez = primers;
	if(nrow(Rez) < 1) {
		Rez$GC  = numeric(0);
		Rez$GCp = numeric(0);
		return(Rez);
	}
	tmp = Rez[1,]; nS = tmp$nS;
	nn  = x[seq(nS, length.out = tmp$Len)];
	Rez$GC = -1;
	tmpGC  = sum(nn %in% c("G", "C"));
	Rez$GC[1] = tmpGC;
	if(nrow(Rez) == 1) {
		Rez$GCp = Rez$GC / Rez$Len;
		return(Rez);
	}
	# Multiple Primer-Seq:
	for(nr in seq(2, nrow(Rez))) {
		nS2 = Rez$nS[nr];
		if(nS == nS2) {
			if(x[nS + Rez$Len[nr] - 1] %in% c("G", "C")) tmpGC = tmpGC + 1;
		} else {
			nS = nS2;
			nn = x[seq(nS, length.out = Rez$Len[nr])];
			tmpGC = sum(nn %in% c("G", "C"));
		}
		Rez$GC[nr] = tmpGC;
	}
	Rez$GCp = Rez$GC / Rez$Len;
	return(Rez);
}

as.df.primer = function(x) {
	p = lapply(x, function(x) {
		x = data.frame(x);
		x$Len = seq(x$Len[1], length.out = nrow(x));
		return(x);
	});
	p = do.call(rbind, p);
	class(p) = c("primer.df", class(p));
	return(p);
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
find.simTm = function(nn.seq, Tm, iter = 2000, tol = 0.1) {
	LEN = length(nn.seq);
	for(i in seq(iter)) {
		nn = sample(nn.seq, LEN);
		Tv = Tm.nnSeq(nn);
		if(abs(Tv - Tm) <= tol) return(nn);
	}
	return(NA);
}

### Note:
# Examples moved to file:
# PCR.Examples.R;

