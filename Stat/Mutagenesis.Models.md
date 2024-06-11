
# Modelling Genetic Drift

Ideas how to model the genetic drift in a population.


## Population

- Diploid organism: 23 pairs of chromosomes;
- Females are XX, males are XY;
- assuming equal size of all chromosomes and equal risk of mutation;
- chromosomal recombination can happen: chromosomes are split on average in 5 segments which can be swapped;
- Note: different variants of models can be simulated;

## Fitness

- individuals without any mutation are considered fir;
- lifespan of fit individuals: 40 years;
- fitness decreases with number of mutations: e.g. decreased lifespan;
- TODO: various subtypes of fitness functions;

## Habitat

- population is spread over a 2D area of a certain (large) size;
- population is considered mainly sedentary with a constant density;
- TODO: find/describe mechanism which maintains constant density: e.g. spreading outwards?

## Breeding

- individuals are fertile between 15 years up to 40 years;
- individuals mate each year once;
- the mating partner is selected randomly each year;
- individuals can find a mating partner up to a distance d from their location;
- spontaneous abortion / early offspring death increases exponentially above a certain threshold of mutations in the offspring;
- TODO: clarify additional details of the mating process (e.g. if a discrete modelling approach is used);

## Questions

1. What is the dynamics of mutations in this population?
2. Do spatially distant populations drift away?
3. Beyond what distance and over what time-scales are speciation events likely to happen?
4. Define a meaningful measure of a speciation event: proportion of offspring death (e.g. 50% death rate)?


## References

- various materials which may be interesting or useful;

1. Matthias Birkner: Coalescence times of ancestral lineages in two-dimensional logistic branching random walks
  > https://www.youtube.com/watch?v=VQMt9tsq-9g
