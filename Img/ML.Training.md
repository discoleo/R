
# Image Processing

## Sample Size

Many AI specialists want to train their machine-learning infrastructures using sets of images. They obviously want to know how many images to include in the training set, i.e. they want to perform a sample size calculation.

### Variance of "Population" of Images

In order to compute the sample size, one needs to estimate the variance of the "population" of images.

Several approaches are possible. More general approaches may yield humongous variances, destroying all hope to ever train an AI.

**Wasserstein Distance:**

A more useful approach is based on the Wasserstein distance. The images can be focused on a specific topic and the barycenter can be computed for the set of images.

The Wasserstein distances compared to the barycenter image can be used as a measure of variance, providing a rough estimate of the "population" variance.

