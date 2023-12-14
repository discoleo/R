
# Uniprot: Cleanup

**Consistency checks on the human proteome**

### A0A3B3ITW3 vs P28288

- Genome: ENST00000647998.2;
- potential isoform: predicted;
- same global sequence;
- 12 "polymorphisms"/splicing isoform:
> S116G; {N118S, K119R, A120K, L121D},
> {K124R, P125Y, F126L, F127L, K128N},
> T133M, M135L;

Entry	Chromosome	Entry Name	Protein names	Gene Names	Length	GeneID	Reviewed\
A0A3B3ITW3	1	A0A3B3ITW3_HUMAN	ATP binding cassette subfamily D member 3	ABCD3	659		unreviewed\
P28288	1	ABCD3_HUMAN	ATP-binding cassette sub-family D member 3 (EC 3.1.2.-) (EC 7.6.2.-) (70 kDa peroxisomal membrane protein) (PMP70)	ABCD3 PMP70 PXMP1	659	5825;	reviewed


### C9K0E9 vs P22760

- Genome: ENST00000488869.1, mRNA: TSL2;
- Fragment of 204 AA vs 399 AA;
- + L202V, L203C, D204S: possible "cloning" (TSL2) artefact?

C9K0E9	3	C9K0E9_HUMAN	Arylacetamide deacetylase	AADAC	204		unreviewed\
P22760	3	AAAD_HUMAN	Arylacetamide deacetylase (EC 3.1.1.3)	AADAC DAC	399	13;	reviewed


### A0A3B3IRV8 vs P78363

- 3 sub-fragments joined together;
- Fragment 1: 1..720 AA vs 1..720 (length 720 AA);
- Fragment 2: 721..898 vs 795..972 (length 178 AA);
- Fragment 3: 899..979 = Unknown 79 AA;

- possible alternative splicing?

A0A3B3IRV8	1	A0A3B3IRV8_HUMAN	ATP binding cassette subfamily A member 4	ABCA4	979		unreviewed\
P78363	1	ABCA4_HUMAN	Retinal-specific phospholipid-transporting ATPase ABCA4 (EC 7.6.2.1) (ATP-binding cassette sub-family A member 4) (RIM ABC transporter) (RIM proteinv) (RmP) (Retinal-specific ATP-binding cassette transporter) (Stargardt disease protein)	ABCA4 ABCR	2273	24;	reviewed

