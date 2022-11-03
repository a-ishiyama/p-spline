
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
  cat('Order', m$b.order,'p-spline with order',m$p.order, "penalty", "\n")
  cat("Effective degrees of freedom:",m$edf,"   Coefficients:",
      m$number.of.coefficients, "\n")
  cat("residual std dev:",m$sig2^0.5,"   r-squared:",m$r2, "   GCV:", m$gcv)
  
  print_list<-list(m$gcv, m$edf, m$r2)
  invisible(print_list)
}



#Q3

predict.pspline <- function(m,x,se=TRUE){
  x<-x[order(x)] #Order the new x data in an ascending order.
  Xp<- splines::splineDesign(knots=m$knots,x,ord=m$b.order+1,outer.ok=TRUE)
  prediction <- Xp%*%m$coef
  # Save fitted values obtained by applying the estimated coefficients from 
  # the training data and the new Xp observation matrix.
  
  if (se==FALSE){
    list(prediction=prediction)
  }else{
    sigma_matrix <- diag(m$sig2,m$number.of.coefficients)
    cholXp <- chol(crossprod(Xp)+(m$lsp*crossprod(m$D)))
    covariance_matrix <-backsolve(cholXp, forwardsolve(t(cholXp),sigma_matrix))
    se <- rowSums(Xp*(Xp%*%covariance_matrix))^0.5
    list(fit=prediction,se=se)
  }
}




#Q4
plot.pspline <- function(m){
  #
  
  x_4 <- seq(m$x[1],m$x[length(m$x)], length=length(m$x))
  
  upper <- lower <-c()
  # Create vectors containing observation numbers of 0s to later save the upper
  # and lower bound of the confidence interval for the fitted value.

  P <- predict.pspline(m, x_4)
  

  upper <- P$fit + 1.96* P$se
  # Save the upper confidence interval for the estimated fitted value.
  lower <- P$fit - 1.96* P$se
  # Save the lower confidence interval for the estimated fitted value.
    
  
  plot(m$x,m$y,main="Original data with the estimated smooth function", 
       xlab="x", ylab="y")
  # Plot a scatter plot of the original y data against x.
  lines(x_4,P$fit, col="red")
  # Draw the estimated smooth function inside the above plot in red.
  lines(x_4,upper, lty=3, col="blue")
  # Draw the upper confidence interval in a blue dotted line in the above plot.
  lines(x_4,lower, lty=3, col="blue")
  # Draw the lower confidence interval in a blue dotted line in the above plot.
  legend(30, -95, legend=c("estimated smooth function",
                           "95% credible intervals"),col=c("red", "blue"), 
         lty=1:2)
  # Set legend to increase the visibility.
  
  plot(m$fitted, m$residual, main="Model residuals against fitted values", 
       xlab="fitted value", ylab="model residuals")
  # Plot the residuals against the estimated fitted value from pspline.
  abline(h=0) 
  # Set a line at y=0 to visibly assess the assumption of residual normality.
  
  qqnorm(m$residual) # Create the Q-Q plot using the residual from pspline.
  qqline(m$residual) # Set a reference line in the above Q-Q plot.
  
  plot_list<-list(cbind(ll=lower,ul=upper,x=x_4))
  invisible(plot_list)
  # Invisibly return a list containing the lower and upper confidence interval,
  # and the x observations.
}