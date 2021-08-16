# fast-pca
Principal Component Analysis implemented in C++

Principal Component Analysis (PCA) is a popular dimensionality reduction algorithm. PCA identifies the hyperplane that lies closest to the data and then projects the data onto it. The hyperplane is given by the eigenvectors of the covariance matrix of the data. The corresponding eigenvalues are proportional to the explained variance.

## Eigendecomposition

Let $\textbf{A}$ be a square n x n matrix with n linearly independent eigenvectors. Then $\textbf{A}$ can be factorised as:

$$\textbf{Av} = \lambda \textbf{v}$$
$$\textbf{AQ} = \textbf{Q} \Lambda$$
$$\textbf{A} = \textbf{Q} \Lambda \textbf{Q}^{-1}$$

This repo contains code which implements the power-deflate method to find eigenvectors and eigenvalues of a data set. 

### The Power Method
The power method is a numerical technique to estimate the dominant eigenpair of a matrix. An eigenvalue $\lambda_k$ is dominant if $|\lambda _ k | \ge |\lambda _ j| \ \ \ \ \ \forall j$. We refer to such $(\textbf{s}_k , \lambda_k)$ as dominant eigenpairs.

Pseudocode:
```matlab
In: s % initial eigenvector iterate

s <- s/||s||    % normalise
while not converged do
    t = As
    |lambda| = ||t|| % eigenvalue estimate
    s = t/|lambda|   % normalised eigenvalue estimate
end while
if (As)/s < 0 then
    lambda_hat <- -lambda
else
    lambda_hat <- lambda

Out: (s, lambda_hat)   % converged eigenpair estimate
```

### The Deflate Method
Deflation removes the largest eigenvector from the matrix.

$$\textbf{C} = \textbf{A} - \lambda_1x_1x_1^T$$

Where $x_1$ is the normalised eigenvector corresponding to largest eigenvalue $\lambda_1$.

Combining this with the power method, we can numerically compute successive eigenpairs of a matrix.

These methods have been implemented in C++ code in the file `myEigenFunctions.cpp`

## Principal Components
PCA assumes that each column of the feature matrix has zero mean. As such, the covariance matrix must be computed only after *centering* the feature matrix; ie ensuring that each column has a mean of 0.

$$\text{cov} = \frac{1}{n-1} \textbf{X}_c^T \textbf{X}_c $$

Where $ \textbf{X}_c $ is the centered feature matrix.


## The Data

The data set used in this example is a set of spectra which were produced from a Fourier-Transform Infrared Spectrometer (FTIR). The reflectance spectra belong to two different types of coffee beans (Arabica and Robusta). The objective of this example is to demonstrate the utility of PCA to aid in the characterisation and classification of the coffee beans. The raw data is contained in the file `DS19hH2dk0FTIRSpectrainstantcoffee.csv`. 