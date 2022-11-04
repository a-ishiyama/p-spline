# Name: Akira Ishiyama
# ID: s2445245

#   The following codes produce a class and related useful method function for
# estimating a smoothing function with basis expansion and penalties.
#   The class is named 'pspline' which 


pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100){
  #   The following set of codes computes a class 'pspline' which estimates the 
  # estimated best fit smooth line for a set of given data x and y.
  #   hh.

  lsp <- seq(logsp[1],logsp[2],length=ngrid)
  # Create a vector with numbers ranging from the lower end to the
  # upper end of logsp. This vector should contain ngrid length of number
  # equally spaced. This is later used to compute the smoothing parameter.
  edf <- gcv <- c()
  # Save two empty vectors to later insert the obtained effective degree of
  # freedom and generalised cross validation criterion score.
  
  dk <- diff(range(x))/(k-bord)
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1) 
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  D <- diff(diag(k),differences=pord)
  
  idenk <- diag(k) # Create an identity matrix with k number of rows.
  crossD <- crossprod(D) # Cross product of D matrix for later use.
  p<-ncol(X);n <-nrow(X) 
  # Save the column count of the X matrix as p, as the number of coefficients to
  # be estimated.
  # Also, save the row count of the X as n, as the number of observation.
  QR<-qr(X) # QR decomposition of matrix X.
  R<-qr.R(QR) # Save R as the R component from the QR composition.
  
  AR <-forwardsolve(t(R),crossD) 
  # Forwardsolve the transpose of R from the QR decomposition of the X matrix
  # with the cross product of D. This step produces the product needed for the
  # forwardsolve in the next step.
  A <-forwardsolve(t(R), t(AR)) 
  # Forwardsolve the transpose of R with the previously derived product AR.
  # These steps computes the product A, needed for the eigen decomposition. 
  ep <- eigen(A) # Eigen decomposition of the product A.
  evl <-diag(ep$values) # Save the eigenvalue from the eigen decomposition.
  evc <-ep$vectors # Save the eigenvector from the eigen decomposition.
  
  for (i in 1:ngrid){
    # Loop to find the optimal smoothing parameter based on generalised cross 
    # validation criterion (gcv).
    core <- (idenk + exp(lsp[i])*evl)
    # Core here refer to value obtained by the addition of identity matrix with
    # the multiplication of smoothing parameter and the previously calculated
    # eigenvalue. The smoothing parameter is expressed in the exponential form
    # as 
    edf[i]<-sum(1/diag(core))
    # Effective degree of freedom is used instead of conventional degree of
    # freedom as the estimated smooth function can vary from being a very 
    # wiggly to a straight line depending on the smoothing parameter.
    QRnewcore<-qr(t(R)%*%evc%*%core%*%t(evc)%*%R)
    cholnewcore <-qr.R(QRnewcore)
    beta <- backsolve(cholnewcore, qr.qty(QRnewcore,t(X)%*%y)[1:p])
    y_hat <-X%*%beta
    # Fitted value is obtained by multiplying the X matrix with the 
    # coefficients.
    sigma <- sum((y-y_hat)^2)/(n-edf[i])
    # Sigma here saves the variance of the estimation residual for each loop
    # of using different smoothing parameter.
    gcv[i]<-sigma/(n-edf[i])
    # Assign the generalised cross validation criterion score for each smoothing 
    # parameter using the definition of gcv, 
    # variance/(number of observation - edf)
  }
  
  i.opt<-which(gcv==min(gcv))
  # Search for the lowest generalised cross validation criterion score.
  core <- (idenk + exp(lsp[i.opt])*evl)
  # The following steps from "core" to "cholnewcore" essentially repeat the
  # same steps as when searching for the optimal smoothing parameter. 
  # This set of step is repeated to compute the optimal coefficients.
  QRnewcore<-qr(t(R)%*%evc%*%core%*%t(evc)%*%R)
  cholnewcore <-qr.R(QRnewcore)
  beta <- backsolve(cholnewcore, qr.qty(QRnewcore,t(X)%*%y)[1:p])
  # Estimate the optimal coefficient.
  y_hat <-X%*%beta
  # Save the predicted value by multiplying the X data with the corresponding
  # estimated coefficients.
  sigma <- sum((y-y_hat)^2)/(n-edf[i.opt])
  # This derives the variance (sigma) of the estimation residual.
  r2 <- 1 - ((n-1)*sigma)/sum((y-mean(y))^2)
  # Derive the r-squared from the estimation and save it for print function.
  residual <- y-y_hat # Safe the residuals.
  
  best_fit_list <- list(coef=beta,fitted=y_hat,sig2=sigma,edf=edf[i.opt],
                        gcv=gcv[i.opt],knots=knots,number.of.coefficients=k, 
                        r2=r2, b.order=bord, p.order=pord, residual=residual, 
                        x=x, y=y, lsp=exp(lsp[i.opt]), D=D)
  # Save information that is critical in the following construction of the 
  # method functions as 'best_fit_line'.
  class(best_fit_list)="pspline" 
  # Defines the class of what is to be returned from this function as 'pspline'
  best_fit_list 
}



print.pspline<- function(m){
  #   This function defines the method function print() for the class 'pspline'.
  #   It prints the essential information in understanding the estimated model,
  # including the order of B-spline and penalty, effective degree of freedom,
  # number of estimated coefficients, residual standard deviation, r-squared,
  # and the generalised cross validation criterion. 
  #   It also silently returns a list containing gcv, edf, and r2.
  
  cat('Order', m$b.order,'p-spline with order',m$p.order, "penalty", "\n")
  # These three lines prints out the aforementioned information.
  cat("Effective degrees of freedom:",m$edf,"   Coefficients:",
      m$number.of.coefficients, "\n")
  cat("residual std dev:",m$sig2^0.5,"   r-squared:",m$r2, "   GCV:", m$gcv)
  
  print_list<-list(m$gcv, m$edf, m$r2)
  invisible(print_list)
  # Invisibly return a list containing the generalised cross validation 
  # criterion (gcv) of the fitted mode, effective degree of freedom (edf),
  # and the r-squared (r2).
}



predict.pspline <- function(m, x, se = TRUE){
  #   The following defines a method function predict() for the class 'pspline'.
  #   It takes the argument se = TRUE/FALSE depending on the user if he/she 
  # wants to obtain the standard error alongside the predicted value.
  #   The prediction is made using the estimated coefficients from m, belonging 
  # to the class 'pspline', and a new set of x values.
  #   The new x data vector is first transformed into an Xp matrix with the
  # number of observation of x as its number of row, and the number of column as
  # the same as the number of coefficients. This Xp matrix is computed the same
  # way as when the estimated coefficients are first estimated for pspline.
  # knots, B-spline and penalty orders are the same as for the training data.
  #    The derivation of covariance matrix for the construction of standard
  # errors is efficiently computed using Cholesky decomposition.
  
  Xp <- splines::splineDesign(knots = m$knots, x, ord = m$b.order+1,
                              outer.ok = TRUE)
  # This creates the X matrix using x vector essentially the same method as
  # used when we first estimated the coefficients etc, in pspline.
  
  prediction <- Xp%*%m$coef
  # Save fitted values obtained by applying the estimated coefficients from 
  # the training data and the new Xp observation matrix.
  
  if (se == FALSE){
    list(prediction = prediction) 
    # If se argument in the method function is set as FALSE, then s.e. is not 
    # required by the user; hence, returning only the list containing the 
    # predictions.
    
  }else{
    # The following lines are built to estimate the standard error of the
    # prediction when the user state the argument for se as TRUE, which is the
    # default.
    
    sigma_matrix <- diag(m$sig2, m$number.of.coefficients)
    # Create a variance matrix that has the variance of the estimation residuals
    # on the diagonal line. The number of column should be the number of
    # coefficients that have been estimated as we are interested in the 
    # covariance of the coefficients.
    
    cholXp <- chol(crossprod(Xp) + (m$lsp * crossprod(m$D)))
    # The covariance matrix has the formula (X^T*X+lambda*D^T*D)^{-1}*sigma^2.
    # It can be inefficient to compute the product before sigma^2. Hence,
    # first Cholesky decompose (X^T*X+lambda*D^T*D).
    covariance_matrix <-backsolve(cholXp, forwardsolve(t(cholXp),sigma_matrix))
    # Then, forwardsolve the transpose of the product from the Cholesky
    # decomposition with sigma^2. This will provide the product to backsolve
    # with the product of Cholesky decomposition.
    
    se <- rowSums(Xp * (Xp %*% covariance_matrix)) ^0.5
    # Save the standard error. This method is more efficient in computation
    # compared to computing (Xp%*%Covariance_matrix%*%t(XP)).
    
    list(fit = prediction, se = se) 
    # Return the list containing the fitted value and the standard error.
  }
}



plot.pspline <- function(m){
  #   The following creates a method function plot() that produces three plots 
  # for the results to be obtained from the estimation done in class 'pspline'.
  #   The first plot produces a scatter plot of the original y data against 
  # x data, and on top of it draws an estimated smooth function line. 
  # This plot also has dotted blue line that represents the 95% confidence 
  # interval for the smooth line.
  #   The second plot produces a residual plot, plotting residuals against the 
  # fitted value from pspline. It is to check whether the model satisfies the
  # assumption that the residual has mean zero with no obvious distribution
  # pattern.
  #   The last plot produces a Q-Q plot of the residuals. It enable the user to 
  # visibly check whether the model satisfies the residual normality 
  # assumption.
  
  upper <- lower <- c()
  # Create empty vectors to later save the upper and lower bound of the 
  # confidence interval for the fitted value.
  
  x_4 <- seq(min(m$x), max(m$x), by=0.1)
  # Equally spaced 
  
  P <- predict.pspline(m, x_4)
  # Use the previously defined method function of predict() to estimate
  # the fitted value and standard error for the computation of the 95%
  # confidence interval.

  upper <- P$fit + 1.96 * P$se
  # Save the upper confidence interval for the estimated fitted value.
  lower <- P$fit - 1.96 * P$se
  # Save the lower confidence interval for the estimated fitted value.
  
  plot(m$x,m$y,main = "Original data with the estimated smooth function", 
       xlab = "x", ylab = "y")
  # Plot a scatter plot of the original y data against x.
  lines(x_4,P$fit, col = "red")
  # Draw the estimated smooth function inside the above plot in red.
  lines(x_4,upper, lty = 3, col = "blue")
  # Draw the upper confidence interval in a blue dotted line in the above plot.
  lines(x_4,lower, lty = 3, col = "blue")
  # Draw the lower confidence interval in a blue dotted line in the above plot.
  legend(30, -95, legend = c("estimated smooth function",
                           "95% credible intervals"),col = c("red", "blue"), 
         lty = 1:2)
  # Set legend to increase the visibility.
  
  plot(m$fitted, m$residual, main = "Model residuals against fitted values", 
       xlab ="fitted value", ylab = "model residuals")
  # Plot the residuals against the estimated fitted value from pspline.
  abline(h=0) 
  # Set a line at y=0 to visibly assess the assumption of residual having a mean
  # value of zero, and randomly scattered.
  
  qqnorm(m$residual) # Create the Q-Q plot using the residual from pspline.
  qqline(m$residual) # Set a reference line in the above Q-Q plot.
  
  plot_list <- list(ll = lower, ul = upper, x = x_4)
  invisible(plot_list)
  # Invisibly return a list containing the lower (ll) and upper (ul) confidence 
  # interval,and the x observations (x).
}


library(MASS)
x<-mcycle[,1]
y<-mcycle[,2]
demo <-pspline(x,y)
#print,predict,plot