# Smoothing with basis expansions and penalties
Author: Akira Ishiyama

This is a repo for the individual project of the module 'Statistical Programming'.

  The R codes produce a class "pspline" and related useful method
function for estimating smoothing function with basis expansion and penalties.


  When given a set of observation containing explanatory variable x and
explained variable y, the user of this class can set their desired order for
B-spline basis function, to approximate the smooth function using basis
function. The user can also  set the number of basis function to use. Given
that this basis approximation is flexible and able to capture a wide range of
function complexities, there is a possibility that the user may fall on the
trap of over-fitting. To avoid such a danger, the user of this class can also
impose a smooth penalty on the coefficients to encourage the estimated fit to
vary smoothly, instead of capturing random noises and turning into a rather
wiggly fit. This class function will automatically selects the best smoothing
parameter (the penalty), by searching for the optimal value based on the
generalised cross validation criterion. Some matrix decomposition methods are
used along to achieve an efficient computation.


  The class pspline comes with three method function, the first being print().
This method function prints some of the crucial information to be obtained
by using the class pspline, such as the residual standard deviation and
r-squared. This helps the user to understand the output in a clean manner.


  The predict() method function enables the user the predict a new set of
fitted values using the estimated coefficients from pspline and a new set of
x observation data. Depending on the users' need, this method function can
also return the standard errors of the fitted values.


  The plot() method function enables the user to plot the original x, y data
with the fitted smooth line and its 95% confidence interval. It also allows
the user to check the assumption made to validate the estimation, namely the
zero mean of the residual and its normality.

