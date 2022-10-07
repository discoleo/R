########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Heterogeneous Symmetric
### Old History
###
### draft v.0.1a

### Old History
# - the historic details;


###############
### History ###
###############

### [branch v.0.3]
#
### v.0.3e & v.0.3g: [cleanup]
# - moved Old History to file:
#   Poly.System.Hetero.Symmetric.S2.History.R;
# - moved section with Multiple Simple Leading terms to file:
#   Poly.System.Hetero.Symmetric.S2.L2.R;
# - moved Multiple/Mixed Leading to file:
#   Poly.System.Hetero.Symmetric.S2.LMixed.R;
# - moved Lnn (Mixed/Trivial Leading) to file:
#   Poly.System.Hetero.Symmetric.S2.Lnn.R;
# - moved section with Mixed Leading Term to new file:
#   Poly.System.Hetero.Symmetric.S2.Leading.R;
# - more refactoring & cleanup;
### v.0.3f:
# - solved: x^4 + b3*y^3 + b2*y^2 + b1*y = R;
### v.0.3e: [refactoring]
### v.0.3c - v.0.3c-clean2:
# - more cleanup;
# - more examples & more A1 extensions;
# - more classic polynomials; [v.0.3c-clean2]
### v.0.3b-ex:
# - more/various examples:
#   5 + 3*x^2 + 6*x^3 + x^6 = 0;
### v.0.3b:
# - Extensions of type A1:
#   x^3 + b1*y + b2*(x+y) = R;
#   Example polynomial: -3 + 3*x^2 - 2*x^3 + x^6 = 0;
### v.0.3a-clean1 - v.0.3a-clean5:
# - started to move derivations to file:
#   Poly.System.Hetero.Symmetric.Derivation.R;
# - more classical polynomials;
### v.0.3a-pre:
# - more work on classic polynomials;
# - TODO: cleanup;


### [branch v.0.1]
### draft v.0.1n:
# - solved:
#   -- x^3 + a*x^3*y = R; (TODO: P3 => P6)
#   -- x^3 + a*x*y^3 = R; (TODO: P5 => P10)
### draft v.0.1m - v.0.1.m-ext:
# - various extensions:
#   -- x^4*y^3 Series;
#   -- x^3*y Series [v.0.1.m-ext];
# - more work on older P6 polynomials;
### draft v.0.1l:
# - solved/extension:
#   a1*x^3 + a2*y^3 + b2*x*y + b1*x = R; (DONE: P3 => P6)
# - improved formatting;
### draft v.0.1k & v.0.1k-x:
# - solved: x^5 + b*y = R;
# - worked out various older issues;
# - P6 for x^3 + b2*x*y + b1*x [DONE in v.0.1k-x];
# - some extensions:
#   O4.1b.) x^4 + b2*x*y + b1*y = R; (DONE: P6 => P12)
#   O4.1c.) x^4 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P6 => P12)
### draft v.0.1j:
# - added a1*x^3*y + a2*x*y^3 + b*x = R;
### draft v.0.1i:
# - added x^3 + b3*x*y + b2*y^2 + b1*y = R;
# - initial exploration of: x^5 + b*(x+y) = 0 variants;
# - more classical polynomials computed;
### draft v.0.1h - v.0.1h-x:
# - initial work on:
#   x^3 + b3*x^2*y + b2*x*y^2 + b1*x = R;
# - improved formatting + some fixes + completion of some cases;
### draft v.0.1f - v.0.1g:
# - added various Mixt-High-Power variants:
#   x^4*y + b*x = R;
#   x^2*y^2 + ... variants; [v.0.1g]
### draft v.0.1e:
# - added shifted version: (x-s)^4 + b*y = R;
### draft v.0.1d:
# - added variant with 2 high-power terms:
#  -- variant 1: a1*x^3 + a2*y^3 + b*x;
#  -- variant 2: a1*x^3 + a2*y^3 + b*x*y;
### draft v.0.1c:
# - added x^3 + b1*x*y + b2*x = R;
# - added x^3 + b1*x*y + b2*y = R;
# - added x^3 + b1*(x*y)^2 = R;
# - TODO:
#  -- shifted versions;
#  -- parametric polynomials [DONE: x & y variants];
### draft v.0.1b - v.0.1b-x:
# - added a basic xy-type: x^3 + x*y = R;
# - added also the shift (v.0.1b-sh);
# - TODO: parametric classic polynomial [DONE: in v.0.1b-x];
### draft v.0.1a-shift:
# - derivation of the classical polynomial for shifted root;
# - more interesting polynomials are generated,
#   when the shifted root is shifted back;
#   [in general not identical to non-shifted root polynomials]
### draft v.0.1a:
# - initial version: basic heterogenous systems;


##################
### S3 Systems ###
##################

### [branch v.0.2] P3 Systems:
# - moved to new file:
#   Poly.System.Hetero.Symmetric.S3.R;
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
# - the simple cases are less rewarding; [but major work since then!]

