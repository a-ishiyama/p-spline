
library(MASS)
x<-mcycle[,1]
y<-mcycle[,2]


#Q1

pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100){
  x<-x[order(x)] # Order x from the data in an ascending order.
  y<-y[order(x)] # Order y from the same data in the SAME order as x.

  lsp <- seq(logsp[1],logsp[2],length=ngrid)
  # 
  edf <- gcv <- c()
  # Save two empty vectors to later insert the obtained edf and gcv.
  
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
  
  AR <-forwardsolve(t(R),crossD) # 
  A <-forwardsolve(t(R), t(AR)) #
  ep <- eigen(A) # Eigen decomposition of the product A.
  evl <-diag(ep$values) # Save the eigenvalue from the eigen decomposition.
  evc <-ep$vectors # Save the eigenvector from the eigen decomposition.
  
  for (i in 1:ngrid){
    core <- (idenk + exp(lsp[i])*evl)
    edf[i]<-sum(1/diag(core))
    QRnewcore<-qr(t(R)%*%evc%*%core%*%t(evc)%*%R)
    cholnewcore <-qr.R(QRnewcore)
    beta <- backsolve(cholnewcore, qr.qty(QRnewcore,t(X)%*%y)[1:p])
    y_hat <-X%*%beta
    sigma <- sum((y-y_hat)^2)/(n-edf[i])
    gcv[i]<-sigma/(n-edf[i])
  }
  
  i.opt<-which(gcv==min(gcv))
  core <- (idenk + exp(lsp[i.opt])*evl)
  QRnewcore<-qr(t(R)%*%evc%*%core%*%t(evc)%*%R)
  cholnewcore <-qr.R(QRnewcore)
  beta <- backsolve(cholnewcore, qr.qty(QRnewcore,t(X)%*%y)[1:p])
  y_hat <-X%*%beta
  sigma <- sum((y-y_hat)^2)/(n-edf[i.opt])
  r2 <- 1 - ((n-1)*sigma)/sum((y-mean(y))^2)
  residual <- y-y_hat
  best_fit_list <- list(coef=beta,fitted=y_hat,sig2=sigma,edf=edf[i.opt],
                        gcv=gcv[i.opt],knots=knots,number.of.coefficients=k, 
                        r2=r2, b.order=bord, p.order=pord, residual=residual, 
                        x=x, y=y, lsp=exp(lsp[i.opt]), D=D)
  class(best_fit_list)="pspline"
  best_fit_list
}

demo <-pspline(x,y)


#Q2

print.pspline<- function(m){
  #   This function defines the method function print() for the class 'pspline'.
  #   
  
  cat('Order', m$b.order,'p-spline with order',m$p.order, "penalty", "\n")
  cat("Effective degrees of freedom:",m$edf,"   Coefficients:",
      m$number.of.coefficients, "\n")
  cat("residual std dev:",m$sig2^0.5,"   r-squared:",m$r2, "   GCV:", m$gcv)
  
  print_list<-list(m$gcv, m$edf, m$r2)
  invisible(print_list)
}



#Q3

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
  
  x <- x[order(x)] #Order the new x data in an ascending order.
  
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




#Q4
plot.pspline <- function(m){
  #   The following creates a method function that produces three plots for 
  # the results to be obtained from the class 'pspline'.
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
  
  x_4 <- seq(m$x[1],m$x[length(m$x)], by = 0.1)
  # 
  
  P <- predict.pspline(m, x_4)
  # Use the previously defined method function of predict.pspline to estimate
  # the 

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
  # Invisibly return a list containing the lower and upper confidence interval,
  # and the x observations.
}