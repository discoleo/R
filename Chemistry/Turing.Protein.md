
# The Turing Protein Gedanken-Experiment

We will explore the reliability of structure-prediction algorithms using the Turing Protein Gedanken-Experiment.

## Setup

The primary protein consists of 3 polypeptide "chains": A, B amd C, forming the linear protein A-B-C. The polypeptide B will be cleaved using endopeptidases, leaving only the polypeptides A & C in the final protein.

- Polypeptide chain A: contains 2 cysteins, Cys-1 & Cys-2 (the indexes 1 & 2 are relative to the cysteins);
- Polypeptide chain C: contains 4 cysteins, Cys-3, up to Cys-6;
- the 2 chains A and C will be linked by 2 cysteine disulfide bonds;
- the chain C will also contain an internal cysteine disulfide bond;
- final protein: A=C=|, where "=|" indicates the intra-chain disulfide bond;

The chains A & C are prespecified and fixed. We will explore multiple chains B: the prerequisite is that endopeptidases cleave off the entire chain B.

## Protein Folding

It is assumed that the amino-acids start to fold as the protein synthesis progresses from the N-terminus:
- as hydrophobic amino-acids emerge from the ribosome channel and come into contact with the aqueous environment of the cytoplasm or endoplasmic reticulum, they start to fold and minimize the surface exposed to the aqueous environment;
- the folding progresses sequentially from N-terminus to C-terminus, as the peptide chain gets elongated;
- the protein folding starts in chain A, then continuous over chain B and terminates with chain C and is therefore dependent on the explicit structure of the chain B;

## Disulfide Bonds

Cysteine Cys-1 can bind to any of the 4 cysteins in chain C, while Cys-2 can bind to any of the 3 remaining cysteins:
- therefore, there are 4 * 3 = 12 possible combinations;
- we will assume that the remaining 2 cysteins are always in suitable positions to form an intra-chain disulfide bond (and further stabilize the conformation);
- we will assume that it is possible to construct suitable polypeptide chains B, each of which will drive the formation of one of the 12 variants; Turing would have found a solution for every variant, as he was quite clever;

## Results

Using our Gedanken-Experiment, we have constructed a set of 12 proteins:
- each one consists of the same polypeptide chains A & C;
- these 2 chains are linked using 2 pairs of disulfide bonds;
- each variant has a different set of disulfide bonds, and it is assumed that the secondary structure of the chains A & C may vary across the 12 protein variants;

## Analysis

Because the amino-acid sequence of the A & C chains is the same, the final proteins will have the same primary amino-acid sequence. The differences are due to the different disulfide bonds.

These differences arise due to the different B chains, which drive different foldings of the A-B-C protein. However, the B chains were cleaved away in the final protein and cannot be used by prediction algorithms which rely solely on the final amino-acid sequence.

## Discussion

Algorithms relying solely on the amino-acid sequence of the final protein may fail miserably to predict the structure of a protein. This may be particularly true when large leading peptides are cleaved away from the final protein.

Another scenario where this can have dramatic effects is in proteins where a middle polypeptide sequence is cleaved away and the free terminals remain linked by disulfide bonds. Although situations with 12 real variants may not be common in nature, the expansion of synthetic peptides may generate such examples.

### Robust Algorithms

New algorithms implementing the sequential folding of the protein, starting from the N-terminus, may be more robust than naive algorithms.

The training of models should use only **entire** proteins, which consist of a single linear amino-acid chain with at most a few amino-acids cleaved from the N-terminus.

Algorithms which fail to incorporate the full structure of the N-terminal subsequence (e.g. focus entirely on the local structure) or which naively base the prediction on arbitrary proteins, may fail miserably in the future, especially as a vast amount of new synthetic polypeptides is generated.

### Spin Glasses

Protein structures may resemble the properties of spin-glasses in various ways:
- a given structure corresponds to a local minimum;
- any meta-stable structure has quasi-"memory": it requires a lot of effort/energy to convert to another structure;

**Explanation:**\
We will use the following assumptions to illustrate a simple case:
- the conversion of structure 1 to structure 2 requires the movement of 1 alpha-helix (or of 1 beta-sheet);
- the corresponding alpha-helix (or beta-sheet) contains 10 amino-acids;
- each amino acid has on average 10 "heavy" atoms (we count only the C, N, O, S atoms which are much heavier than H, and the average is very rough);
- the change in conformation will require the movement of 10 * 10 = 100 atoms;
- thermal fluctuations may be insufficient to drive such changes: alpha-emitting radionuclides may posses the necessary kinetic energy, but most humans do not have such nuclides coordinated to their alpha-helices and beta-sheets;
- proteins can be viewed as behaving like spin-glasses, with their secondary and tertiary structures exhibiting the property of "memory";
- as protein folding progresses from the N-terminus, the existing conformation which is caught into a local minimum will persist due to this spin-glass "memory";

