

### R-Help List:
# Initial Question:
# https://stat.ethz.ch/pipermail/r-help/2023-January/476725.html
# System of Eqs:
# https://stat.ethz.ch/pipermail/r-help/2023-January/476729.html


library(rootSolve)


solve.AcidSpecies = function(x, H, Mg=0.0006, K = 0.12) {
	# with(list(x), { seems NOT to work with multiroot });
	atp = x[1]; adp = x[2]; Pi = x[3];
	pcr = x[4]; cr = x[5]; lactate = x[6];

	### ATP
	hatp <- 10^6.494*H*atp
	hhatp <- 10^3.944*H*hatp
	hhhatp <- 10^1.9*H*hhatp
	hhhhatp  <- 10*H*hhhatp
	mgatp <- 10^4.363*atp*Mg
	mghatp <- 10^2.299*hatp*Mg
	mg2atp <- 10^(-7)*Mg*mgatp; # TODO: check;
	katp <- 10^0.959*atp*K

	### ADP
	hadp <- 10^6.349*adp*H
	hhadp <- 10^3.819*hadp*H
	hhhadp <- 10*H*hhadp
	mgadp <- 10^3.294*Mg*adp
	mghadp <- 10^1.61*Mg*hadp
	mg2adp <- 10*Mg*mgadp
	kadp <- 10^0.82*K*adp

	### Pi
	hpi <- 10^11.616*H*Pi
	hhpi <- 10^6.7*H*hpi
	hhhpi <- 10^1.962*H*hhpi
	mgpi <- 10^3.4*Mg*Pi
	mghpi <- 10^1.946*Mg*hpi
	mghhpi <- 10^1.19*Mg*hhpi
	kpi <- 10^0.6*K*Pi
	khpi <- 10^1.218*K*hpi
	khhpi <- 10^-0.2*K*hhpi

	### P-Crea
	hpcr <- 10^14.3*H*pcr
	hhpcr <- 10^4.5*H*hpcr
	hhhpcr <- 10^2.7*H*hhpcr
	hhhhpcr <- 100*H*hhhpcr
	mghpcr <- 10^1.6*Mg*hpcr
	kpcr <- 10^0.74*K*pcr
	khpcr <- 10^0.31*K*hpcr
	khhpcr <- 10^-0.13*K*hhpcr

	### Crea:
	hcr <- 10^14.3 * H * cr
	hhcr <- 10^2.512 * H * hcr

	### Lactate
	hlactate <- 10^3.66 * H * lactate
	mglactate <- 10^0.93 * Mg * lactate

	### Total:
	tatp <- atp + hatp + hhatp + hhhatp + mgatp + mghatp + mg2atp + katp
	tadp <- adp + hadp + hhadp + hhhadp + mghadp + mgadp + mg2adp + kadp
	tpi  <- Pi + hpi + hhpi + hhhpi + mgpi + mghpi + mghhpi + kpi + khpi + khhpi
	tpcr <- pcr + hpcr + hhpcr + hhhpcr + hhhhpcr + mghpcr + kpcr + khpcr + khhpcr
	tcr  <- cr + hcr + hhcr
	tlactate <- lactate + hlactate + mglactate
	# tmg <- Mg + mgatp + mghatp + mg2atp + mgadp + mghadp + mg2adp + mgpi + 
	#	+ kghpi + mghhpi + mghpcr + mglactate
	# tk <- K + katp + kadp + kpi + khpi + khhpi + kpcr + khpcr + khhpcr

	total = c(
		tatp - 0.008,
		tadp - 0.00001,
		tpi - 0.003,
		tpcr - 0.042,
		tcr - 0.004,
		tlactate - 0.005);
	return(total);
	# })
}

### Conditions

# free K and Mg constrained to be fixed
#
# Mg <- 0.0006
# K <- 0.12


x0 = c(
	atp = 0.008,
	adp = 0.00001,
	Pi = 0.003,
	pcr = 0.042,
	cr = 0.004,
	lactate = 0.005
) / 30;
# [solved] positive value
# x0[1] = 1E-6;

x = multiroot(solve.AcidSpecies, x0, H = 4E-8)

print(x)

