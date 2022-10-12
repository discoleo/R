########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Heterogeneous Symmetric
### Old History
###
### draft v.0.1a


### Hetero-Symmetric
### Polynomial Systems: 3 Variables



####################

###############
### History ###
###############


### draft v.0.4g:
# - robust solutions for S3P3 Simple & S3P2-Asymmetric Sum;
# - various cleanup;
### draft v.0.4f-clean 1 & 2 & 3:
# - more cleaning;
# - S3P3 MixedSideChain: factorized P[12] => P[8]
#   => 3*8 = 24 true solutions; [v.0.4f-TrueSol]
### draft v.0.4e:
# - cleanup: moved to new file
#   Poly.System.Hetero.Symmetric.S3.Derivation.R;
# - work on x^3 + b2*y^2 + b1*z = R;
### draft v.0.4c - v.0.4d:
# - solved Order 3:
#   x^3 + b*y*z = R;
# - TODO: Case x == y, but != z;
# - moved to new file:
#   Poly.System.Hetero.Symmetric.S3.YZ.R; [draft v.0.4d]
### draft v.0.4b - v.0.4b-full:
# - [full] solution for:
#   x^3 + y^2 = R; [partial solution]
#   x^3 + b2*y^2 + b1*y = R; [v.0.4b-full]
### draft v.0.4a:
# - moved to new file:
#   systems with Composite Leading Term: e.g. x^m*y^n;
# - file: Poly.System.Hetero.Symmetric.S3.Leading.R;
### draft v.0.3h - v.0.3h-fix:
# - Order 3: (partially solved + partial bug fix)
#   x^3 + a1*y^3 + a2*z^3 = R;
### draft v.0.3g - v.0.3g-fix:
# - solved (with extensions of type A1):
#   x^2 + a1*y^2 + a2*z^2 = R;
# - more exploration of this system; [v.0.3g-expl]
# - fixed / full equation; [v.0.3g-fix]
### draft v.0.3f:
# - Structural Extension:
#   a2*(x*y*z)^2 + a1*x*y*z + x*y + b1*y = R;
### draft v.0.3e:
# - initial work on: x^3 + b2*y^2 + b1*y = R;
### draft v.0.3d - v.0.3d-simple:
# - solved: x*y + b*y = R;
# - Note: only A1-type extensions have distinct solutions;
# - simplification of the correct solution (for A1-extensions);
### draft v.0.3c - v.0.3c-poly-shift:
# - solved: x^2 + y^2 + b1*y = R;
# - added Extensions of type A1; [v.0.3c-ext]
# - added Classical polynomial: P[6]; [v.0.3c-poly & fixed minor bug]
# - special Case P[6]: b.ext[2] = -3;
#   e.g. 11 + 2*x + 5*x^2 - 2*x^3 + x^6 = 0;
# - special Case with shifted P[6]: b5 = 0;
### draft v.0.3b - v.0.3b-P2ext:
# - solved: x^3 + b*y*z = R;
# - reordering of sections & better comments; [v.0.3b-ord]
# - fix of the wrong roots (in an older Shift-x system); [v.0.3b-fix]
# - extension A1 for the P2 system;
# - TODO: cleanup;
### draft v.0.3a-ext:
# - extensions of type A1 for Ht S3P2;
# - simplification of the base Eq for Ht S3P2;
### draft v.0.3a-pre:
# - moved Difference types to new file:
#   Poly.System.Hetero.Symmetric.S3.Diff.R;
### draft v.0.2b - v.0.2b-ext: [MOVED]
# - first concepts / solved [v.0.2b-sol]:
#   x^2 - y^2 + b*x*y = R;
# - solved extension A1 (Pow 1): + b[2]*(x+y+z); [v.0.2b-ext]
### draft v.0.2a - v.0.2a-poly:
# - solved: Ht[3, 3, 1];
# - simplified solution: from P11 to P8; [v.0.2a-simplify]
# - classic polynomial: P24; [v.0.2a-poly]
# - TODO: find robust solution;
### draft v.0.1d - v.0.1d-poly:
# - minor fix: in Ht[3, 2, 1];
# - classic polynomial (P6) for Ht[3, 2, 1]; [v.0.1d-poly]
### draft v.0.1c-move:
# - moved Hetero-Mixt (v.0.1c) to separate file;
### draft v.0.1c-pre-alpha - v.0.1c-exact:
# - first look & solved + exact/robust solution: [v.0.1c-exact]
#   x*y^2 + y*z^2 + z*x^2 = R1;
#   [moved to file: Poly.System.Hetero.Symmetric.S3.Mixt.R]
### draft v.0.1b - v.0.1b-fix:
# - solved: x[i]^2 + b2*x[j] + b1*x[k];
# - classical Polynomial (P8) (v.0.1b-clP; fixed in v.0.1b-fix);
#  (done: P8 = P2*P6 in v.0.1b-fix)
# - extension: x[i]^2 + s*x + b3*x*y*z + b2*x[j] + b1*x[k]; (v.0.1b-ext)
# - TODO: find/correct (precision) bug vs correct roots;
### draft v.0.1a:
# - moved to new file
#   from Poly.System.Hetero.Symmetric.R;
########################
### former branch v.0.2:
### in Poly.System.Hetero.Symmetric.R
### draft v.0.2.d:
# - more work on x[i]^3 + b*x[i+1] = R;
# - added: x^2*y*z + b*y = R;
# - various formatting improvements;
### draft v.0.2.c:
# - solved: x[i]^2 + s*x[i] + b*x[i+1] = R;
# - TODO: correct various bugs [DONE];
### draft v.0.2b:
# - some exploration of systems with x*y*z terms;
### branch v.0.2a:
# - more work on systems with 3 variables:
#   "proper" implementation of: x[i]^2 + b*x[k];
# - TODO: robust removal of set of wrong solutions;
#   [or avoid getting superfluous solutions ???]
### branch v.0.2a-pre-a:
# - initial work on systems with 3 variables;
# - the simple cases are less rewarding;

