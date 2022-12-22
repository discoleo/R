


### Prions


### 1. PLAAC Algorithm
# PLAAC = Prion-like Amino Acid Composition;
# http://plaac.wi.mit.edu/details
# - based on HMM:
#   https://github.com/whitehead/plaac/blob/master/cli/src/plaac.java


# 2. Toombs, J.A., B. R. McCarty, and E. D. Ross. 2010.
#    Compositional Determinants of Prion Formation in Yeast.
#    Mol Cell Biol. 2010 Jan; 30(1): 319–332.
#    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2798286/
# 3. Alberti, S., R. Halfmann, O. King, A. Kapila, and S. Lindquist. 2009.
#    A systematic survey identifies prions and illuminates sequence features of prionogenic proteins.
#    Cell 137:146-158.


# TODO:

### SUP35 Libs: from Toombs (2010)

### Lib 1:
ppPrionL1 = c(
	"VNIFPYYN", "VTSGSYNT", "ASNIVMNC", "AHTTNMIV",
	"YNCSVNML", "FSIYMPYK", "LLVHSNAI", "WGARQFNI",
	"VTTDILAM", "RRDYLTRF", "STVICGVI", "IHFWPRAP",
	"HSNVSVIH", "TWAPIMVY", "MFQHGIGV", "TRIWNFSG",
	"YHSVEFRI", "TTVNHHFN", "GSLSLQYF", "IFDIANHS",
	"LQPCYCSR", "MLSSNFIH", "SSGPLNFI", "CLSPAECR",
	"QFVARVFR", "LKSVITWN", "SVHVNSTS");

ppNotPrL1 = c(
	"TDPWVPHP", "NPEVPNAN", "THHSHTLP", "YLPFMDTP",
	"PPIVKPRT", "VDDRHMFS", "CKSVCNFD", "GISTRSQE",
	"VSLSKNRL", "LRDPDTCS", "RKATDLFP", "TAYVRHID",
	"DRYKGKPH", "DPNAALVF", "HIHPLFIH", "TLARRDPP",
	"PNASGIHY", "ADSASNAS", "NGPAYPLA", "SVNPALYR",
	"SGVSTAVR", "LNRITLRN", "IVPRNVNC", "NISPFSKD",
	"MTQNPHIF", "LSARPLGH", "LGNPTFHY", "AQDSHPDI",
	"NNPQYLFK", "DERPWCPE", "GPTMNNRD", "THRHNKHR",
	"KGSPSTPT", "EAPSKSAQ", "RPERRSNP", "ICWHTEPY",
	"CIKHINSI", "PVPSSSQP", "GANSAITN", "SHLWRRNR",
	"DSHTGTPR", "STVPPPHH", "VNCARGTA", "QVASQNGR",
	"SSNKFMHT", "GFTKALPG", "ALSSRQWS", "IDKNLMSH",
	"CFLRSYMG", "VALIPKTA", "HNLANHSH", "KMTTNTKH")


### Lib 2:
ppPrionL2 = c(
	"FANHAHWV", "GTTYAPLF", "WNAFSTYS", "HTVHHIYP",
	"LNTFPHSY", "DIMTNNAE", "SQDYSSYD", "CINTGLWL",
	"HLHMSMLS", "DRHYFAGS", "GGPIFNTK", "SFMAVETR",
	"TWDGIGYR", "SPPFETSP", "GVNTHTSY", "SIHMRVSS",
	"HNDRTAFM", "PQNQTWAD", "PDYFFHPT", "HVPSPAHQ",
	"DSDHHFWP", "TSNTIIRA", "DCLGYPGL", "SMHNGTHR",
	"ESILWASQ", "PRLTNHSS", "FWMQRNSC", "SFSYVTFP",
	"CQINWRTA", "GPPFPGQN", "VASWASVG", "YREGDNLW",
	"HTLVFNDR")

ppNotPrL2 = c(
	"SVSDHTNP", "KGRVSGPE", "ATSPVPRH", "YEYSPLQH",
	"TMTDLPYL", "ESILWASQ", "FTRAKSRT", "TTSYHPEL",
	"VAHCRHPL", "SSTLLDPK", "IETHFTLS", "APHGLGPT",
	"RCSDSQGV", "VHHDPVST", "HPIMSSLS", "LGPVHYRN",
	"SMHNGTHR", "DGPTYDWT", "PYKAATRN", "PTYNDPST",
	"LSQSYVQE", "YDSGTPPK", "SQQRFNPT", "HRDNCRTR",
	"PPQAVYPP", "QHASGRDG", "QTRFYGIH", "QTTTAIHA",
	"PHEAVSSC", "RRHYAPSI", "KYMYHANM", "LADSNTPR",
	"SLAAPRDN", "FWIDGSAD", "DRHYFAGS", "IRTHMSSK",
	"ARNMTRYL", "RAYDILPV", "NEDPGTDT", "SRSIRYDN",
	"SQDYSSYD")

