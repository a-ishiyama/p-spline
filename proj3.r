
library(MASS)
x<-mcycle[,1]
y<-mcycle[,2]


#Q1

pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100){
  x<-x[order(x)]
  y<-y[order(x)]
    
  lsp <- seq(logsp[1],logsp[2],length=ngrid)
  edf <- gcv <- rep(0,ngrid)
  
  dk <- diff(range(x))/(k-bord)
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1) 
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  D <- diff(diag(k),differences=pord)
  
  idenk <- diag(k)
  crossD <- crossprod(D)
  p<-ncol(X);n <-nrow(X)
  QR<-qr(X)
  R<-qr.R(QR)
  
  AR <-forwardsolve(t(R),crossD)
  A <-forwardsolve(t(R), t(AR))
  ep <- eigen(A)
  evl <-diag(ep$values)
  evc <-ep$vectors
  
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
  best_fit_list <- list(coef=beta,fitted=y_hat,sig2=sigma,edf=edf[i.opt],gcv=gcv[i.opt],knots=knots,number.of.coefficients=k, r2=r2, b.order=bord, p.order=pord, residual=residual, x=x, y=y, lsp=exp(lsp[i.opt]), D=D)
  class(best_fit_list)="pspline"
  best_fit_list
}

ok <-pspline(x,y)


#Q2

print.pspline<- function(m){
  cat('Order', m$b.order,'p-spline with order',m$p.order, "penalty", "\n")
  cat("Effective degrees of freedom:",m$edf,"   Coefficients:",m$number.of.coefficients, "\n")
  cat("residual std dev:",m$sig2^0.5,"   r-squared:",m$r2, "   GCV:", m$gcv)
  
  print_list<-list(m$gcv, m$edf, m$r2)
  invisible(print_list)
}

print(ok)


#Q3

predict.pspline <- function(m,x,se=TRUE){
  x<-x[order(x)]
  Xp<- splines::splineDesign(knots=m$knots,x,ord=m$b.order+1,outer.ok=TRUE)
  prediction <- Xp%*%m$coef
  
  if (se==FALSE){
    list(prediction=prediction)
  }else{
    sigma_matrix <- diag(m$sig2,m$number.of.coefficients)
    cholXp <- chol(crossprod(Xp)+(m$lsp*crossprod(m$D)))
    covariance_matrix <-backsolve(cholXp, forwardsolve(t(cholXp),sigma_matrix))
    se <- rowSums(Xp*(Xp%*%covariance_matrix))^0.5
    list(prediction=prediction,se=se)
  }
}

predict(ok,x)



#Q4
plot.pspline <- function(m){
  upper <-c(0,length(m$y))
  lower <-c(0,length(m$y))
  for (i in 1:length(m$y)){
    t_value <- qt(0.975, df=length(m$y)-m$number.of.coefficients)
    upper[i] <- m$fitted[i] + t_value* predict.pspline(m, m$x)$se[i]
    lower[i] <- m$fitted[i] - t_value* predict.pspline(m, m$x)$se[i]
  }
  
  plot(m$x,m$y,main="Original data with the estimated smooth function", xlab="x", ylab="y")
  lines(m$x,m$fitted, col="red")
  lines(m$x,upper, lty=3, col="blue")
  lines(m$x,lower, lty=3, col="blue")
  legend(35, -95, legend=c("estimated smooth function", "95% credible intervals"),col=c("red", "blue"), lty=1:2, cex=0.8)
  
  plot(m$fitted, m$residual, main="Model residuals against fitted values", xlab="fitted value", ylab="model residuals")
  abline(h=0)
  
  qqnorm(m$residual)
  qqline(m$residual)
  
  plot_list<-list(cbind(ll=lower,ul=upper,x=m$x))
  invisible(plot_list)
}

plot(ok)


