
# Image Processing

## Sample Size

Many AI specialists want to train their machine-learning infrastructures using sets of images. They obviously want to know how many images to include in the training set, i.e. they want to perform a sample size calculation.


### Variance of "Population" of Images

In order to compute the sample size, one needs to estimate the variance of the "population" of images.

Several approaches are possible. More general approaches may yield humongous variances, destroying all hope to ever train an AI.

#### Wasserstein Distance:

A more useful approach is based on the Wasserstein distance. The images in the set can be focused on a specific topic and the barycenter can be computed for the specific set of images.

The Wasserstein distances computed based on the barycenter image can be used as a measure of variance, providing a rough estimate of the "population" variance.

**Problems:**\
Unfortunately, the Wasserstein distance has some shortcomings as well. It may transport a density into several locations.

A better approach would preserve the connectivity of the components with the addition of 2 new operations (which do not preserve the mass):
- Amplification or attenuation of a certain density;
- Dilation or shrinkage of a certain area (or density);


**Operations:**\
- Translation;
- Skew;
- Dilation or Shrinkage;
- Amplification or Attenuation;

