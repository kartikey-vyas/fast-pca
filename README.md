# fast-pca
Principal Component Analysis implemented in C++

Principal Component Analysis (PCA) is a popular dimensionality reduction algorithm. PCA identifies the hyperplane that lies closest to the data and then projects the data onto it. The hyperplane is given by the eigenvectors of the covariance matrix of the data. The corresponding eigenvalues are proportional to the explained variance.

## Eigendecomposition

Let A be a square n x n matrix with n linearly independent eigenvectors. Then A can be factorised as:

![equation](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BAv%7D%20%3D%20%5Clambda%20%5Ctextbf%7Bv%7D%20%5C%5C%20%5Ctextbf%7BAQ%7D%20%3D%20%5Ctextbf%7BQ%7D%20%5CLambda%20%5C%5C%20%5Ctextbf%7BA%7D%20%3D%20%5Ctextbf%7BQ%7D%20%5CLambda%20%5Ctextbf%7BQ%7D%5E%7B-1%7D%20%5C%5C)

This repo contains code which implements the power-deflate method to find eigenvectors and eigenvalues of a data set. 

### The Power Method
The power method is a numerical technique to estimate the dominant eigenpair of a matrix.

![equation](https://latex.codecogs.com/gif.latex?%5Ctextup%7BAn%20eigenvalue%7D%20%5C%20%5Clambda_k%20%5C%20%5Ctextup%7Bis%20dominant%20if%7D%20%5C%20%7C%5Clambda%20_%20k%20%7C%20%5Cge%20%7C%5Clambda%20_%20j%7C%20%5C%20%5C%20%5C%20%5C%20%5C%20%5Cforall%20j.%20%5C%5C%20%5Ctextup%7BWe%20refer%20to%20such%7D%20%5C%20%28%5Ctextbf%7Bs%7D_k%20%2C%20%5Clambda_k%29%20%5C%20%5Ctextup%7Bas%20dominant%20eigenpairs.%7D)

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

![equation](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BC%7D%20%3D%20%5Ctextbf%7BA%7D%20-%20%5Clambda_1x_1x_1%5ET)

Where x_1 is the normalised eigenvector corresponding to largest eigenvalue.

Combining this with the power method, we can numerically compute successive eigenpairs of a matrix.

These methods have been implemented in C++ code in the file `myEigenFunctions.cpp`

## Principal Components
PCA assumes that each column of the feature matrix has zero mean. As such, the covariance matrix must be computed only after *centering* the feature matrix; ie ensuring that each column has a mean of 0.

![eqn](https://latex.codecogs.com/gif.latex?%5Ctext%7Bcov%7D%20%3D%20%5Cfrac%7B1%7D%7Bn-1%7D%20%5Ctextbf%7BX%7D_c%5ET%20%5Ctextbf%7BX%7D_c)

Where X_c is the centered feature matrix.


## The Data

The data set used in this example is a set of spectra which were produced from a Fourier-Transform Infrared Spectrometer (FTIR). The reflectance spectra belong to two different types of coffee beans (Arabica and Robusta). The objective of this example is to demonstrate the utility of PCA to aid in the characterisation and classification of the coffee beans. The raw data is contained in the file `DS19hH2dk0FTIRSpectrainstantcoffee.csv`. 
