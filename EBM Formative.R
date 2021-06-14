# This is the empirical example for the paper "Quantile Regression under Misspecification" by
# J. Angrist, V. Chernozhukov and I. Fernandez-Val


########################################################################################################
# This program compute uniform and pointwise confidence bands for the intercept and the schooling coeff.
# For the intercept
# 	The categories of education:    each year
# 	Distribution of the covariates: each year
########################################################################################################

############

############## computes the estimate of the sigma components of the variance formula ####################
##############			(robust to misspecification)			     ####################

sigma <- function(data, n, tau, res)
{
  weight <- perwt * (tau - (res <= 0));
  x <- cbind(weight,weight*educ,weight*exper,weight*exper2,weight*black);
  return(t(x) %*% x /n)
}

sigma2 <- function(data, n, tau, res)
{
  weight <- perwt * sqrt(tau * (1 - tau));
  x <- cbind(weight);
  return(t(x) %*% x /n)
}

############## computes the estimate of the sigma components of the variance formula ####################
##############			(under correct specification)			     ####################

sigma0 <- function(data, n, tau, res)
{
  weight <- perwt * sqrt(tau * (1 - tau));
  x <- cbind(weight,weight*educ,weight*exper,weight*exper2,weight*black);
  return(t(x) %*% x /n)
}


############## computes the estimate of the jacobian components of the variance formula ####################


jacobian <- function(data, n, tau, res, alpha)
{
  hn <- qnorm(1-alpha/2)^(2/3)*(1.5*dnorm(qnorm(tau))^2 / (2*qnorm(tau)^2 + 1))^(1/3)*n^(-1/3)
  weight <- sqrt(perwt)*(abs(res) <= hn)
  x <- cbind(weight,weight*educ,weight*exper,weight*exper2,weight*black);
  return(t(x) %*% x / (2*n*hn))
}

jacobian2 <- function(data, n, tau, res, alpha)
{
  hn <- qnorm(1-alpha/2)^(2/3)*(1.5*dnorm(qnorm(tau))^2 / (2*qnorm(tau)^2 + 1))^(1/3)*n^(-1/3)
  weight <- sqrt(perwt)*(abs(res) <= hn)
  x <- cbind(weight);
  return(t(x) %*% x / (2*n*hn))
}

############## obtains the distribution of K by subsampling ##############################################

subsamplek <- function(data=data, n=n, b=b, B=B, formula, R=R, V, tau, coeffs)
{
  k<-NULL;
  RVR <- (t(R) %*% V %*% R / b)^(-1/2);
  for (s in 1:B)
  {
    sing = 0;
    while (sing == 0)
    {
      sdata = data[sample(n,b,replace=T,prob=perwt), ];
      x = cbind(sdata$perw*sdata$educ,sdata$perw*sdata$exper,sdata$perw*sdata$exper2,sdata$perw*sdata$black,sdata$perw);
      sing = det(t(x) %*% x);
    }
    sqr <- rq(formula, tau, sdata, weights = perwt, method="br");
    k <- cbind(k, abs(RVR %*% t(R) %*% (coeffs - coefficients(sqr))))
  }
  return(k);
}




############# obtains estimates of the coefficients and of the standard errors ###############################
#############			(robust to misspecification)		       ###############################

table.rq.res<-
  function (formula, taus, method = "fn", data=data,
            R=R, alpha=0.05, n, ...)
  {
    m <- length(taus)
    tab <- NULL
    setab<- NULL
    for (i in 1:m) {
      fit <- rq(formula, taus[i], method = "fn", data=data, weights=perwt)
      coeff <- t(R) %*% coefficients(fit)
      tab <- rbind(tab, coeff)
      sigmatau <- sigma(data, n, taus[i], residuals(fit));
      jacobtau <- jacobian(data, n, taus[i], residuals(fit), alpha);
      V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau) / n;
      secoeff <- sqrt(t(R) %*% V %*% R);
      setab <- rbind(setab, secoeff)
    }
    results <- list("tab"=tab, "setab"=setab)
    results
    
  }


table.rq.res2<-
  function (formula, taus, method = "fn", data=data,
            R=R, alpha=0.05, n, ...)
  {
    m <- length(taus)
    tab <- NULL
    setab<- NULL
    for (i in 1:m) {
      fit <- rq(formula, taus[i], method = "fn", data=data, weights=perwt)
      coeff <- t(R) %*% coefficients(fit)
      tab <- rbind(tab, coeff)
      sigmatau <- sigma2(data, n, taus[i], residuals(fit));
      jacobtau <- jacobian2(data, n, taus[i], residuals(fit), alpha);
      V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau) / n;
      secoeff <- sqrt(t(R) %*% V %*% R);
      setab <- rbind(setab, secoeff)
    }
    results <- list("tab"=tab, "setab"=setab)
    return(results)
    
  }


############# obtains estimates of the coefficients and of the standard errors ###############################
#############			(under correct specification)		       ###############################

table0.rq.res<-
  function (formula, taus, method = "fn", data=data,
            R=R, alpha=0.05, n, ...)
  {
    m <- length(taus)
    tab <- NULL
    setab<- NULL
    for (i in 1:m) {
      fit <- rq(formula, taus[i], method = "fn", data=data, weights=perwt)
      coeff <- t(R) %*% coefficients(fit)
      tab <- rbind(tab, coeff)
      sigmatau <- sigma0(data, n, taus[i], residuals(fit));
      jacobtau <- jacobian(data, n, taus[i], residuals(fit), alpha);
      V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau) / n;
      secoeff <- sqrt(t(R) %*% V %*% R);
      setab <- rbind(setab, secoeff)
    }
    results <- list("tab"=tab, "setab"=setab)
    return(results)    
  }



##########################################################################################################################
#######     				I - QUANTILE PROCESS FOR THE INTERCEPT
##########################################################################################################################





##########################################################################################################################
#######     A - SCHOOLING VARIABLE FOR EACH YEAR AND COVARIATES EVALUATED AT EACH YEAR MEAN
##########################################################################################################################


# Part II:  Implementation (Census 80)
setwd ("C:\\Users\\Leandro\\Documents\\BRISTOL\\Modules\\TB2\\EBM\\Computer Classes\\Formative assessment")
library(foreign);
library(quantreg);
memory.limit(size=500000000);

set.seed(8);	

data<-read.dta("census80.dta");
attach(data);

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
alpha<- 0.05;

R<- apply(rbind(1,educ,exper,exper2,black),1,FUN = weighted.mean,w=perwt);

taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20;

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);




detach(data);
rm(data);

##########Saving variables for joint graph ################################################################################################################

res.to.plot.80 <- res.to.plot;
Kalpha.80 <- Kalpha;
olscoeffs.80 <- olscoeffs; 


##########################################################################################################################

# Part II:  Implementation (Census 90)

library(foreign);
library(quantreg);
memory.limit(size=500000000);


set.seed(8);	


data<-read.dta("census90.dta");
attach(data);

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
alpha<- 0.05;

R<- apply(rbind(1,educ,exper,exper2,black),1,FUN = weighted.mean,w=perwt);

taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20;

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);



detach(data);
rm(data);

##########Saving variables for joint graph ################################################################################################################

res.to.plot.90 <- res.to.plot;
Kalpha.90 <- Kalpha;
olscoeffs.90 <- olscoeffs; 


##########################################################################################################################

# Part II:  Implementation (Census 00)

library(foreign);
library(quantreg);
memory.limit(size=500000000);


set.seed(8);	

data<-read.dta("census00.dta");
attach(data);

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
alpha<- 0.05;

R<- apply(rbind(1,educ,exper,exper2,black),1,FUN = weighted.mean,w=perwt);

taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20;

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);


detach(data);
rm(data);

##########Saving variables for joint graph ################################################################################################################

res.to.plot.00 <- res.to.plot;
Kalpha.00 <- Kalpha;
olscoeffs.00 <- olscoeffs; 

##########################################################################################################################

# Part II:  Implementation (Census10)

library(foreign);
library(quantreg);
memory.limit(size=500000000);


set.seed(8);	

data<-read.csv("2010male.csv");
attach(data);
data<- within(data, rm("X"))

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
alpha<- 0.05;

R<- apply(rbind(1,educ,exper,exper2,black),1,FUN = weighted.mean,w=perwt);

taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20;

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);


detach(data);
rm(data);

##########Saving variables for joint graph ################################################################################################################

res.to.plot.10 <- res.to.plot;
Kalpha.10 <- Kalpha;
olscoeffs.10 <- olscoeffs; 

##########################################################################################################################

# Part II:  Implementation (Census19)

library(foreign);
library(quantreg);
memory.limit(size=500000000);


set.seed(8);	

data<-read.csv("2019male.csv");
attach(data);
data<- within(data, rm("X"))

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
alpha<- 0.05;

R<- apply(rbind(1,educ,exper,exper2,black),1,FUN = weighted.mean,w=perwt);

taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20;

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);


detach(data);
rm(data);

##########Saving variables for joint graph ################################################################################################################

res.to.plot.19 <- res.to.plot;
Kalpha.19 <- Kalpha;
olscoeffs.19 <- olscoeffs; 



# Part II:  Implementation WOMEN (Census10)

library(foreign);
library(quantreg);
memory.limit(size=500000000);


set.seed(8);	

data<-read.csv("2010female.csv");
attach(data);
data<- within(data, rm("X"))

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
alpha<- 0.05;

R<- apply(rbind(1,educ,exper,exper2,black),1,FUN = weighted.mean,w=perwt);

taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20;

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);


detach(data);
rm(data);

##########Saving variables for joint graph ################################################################################################################

res.to.plot.fe10 <- res.to.plot;
Kalpha.fe10 <- Kalpha;
olscoeffs.fe10 <- olscoeffs; 





######################## WOMEN CENSUS 2019
library(foreign);
library(quantreg);
memory.limit(size=500000000);


set.seed(8);	

data<-read.csv("2019female.csv");
attach(data);
data<- within(data, rm("X"))

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
alpha<- 0.05;

R<- apply(rbind(1,educ,exper,exper2,black),1,FUN = weighted.mean,w=perwt);

taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20;

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);


detach(data);
rm(data);

##########Saving variables for joint graph ################################################################################################################

res.to.plot.fe19 <- res.to.plot;
Kalpha.fe19 <- Kalpha;
olscoeffs.fe19 <- olscoeffs; 






##########################################################################################################################
#######     			II - QUANTILE PROCESS FOR THE COEFFICIENT OF SCHOOLING
##########################################################################################################################



##########################################################################################################################

# Part II:  Implementation (Census 80)

library(foreign);
library(quantreg);
memory.limit(size=500000000);

set.seed(8);	

data<-read.dta("census80.dta");
attach(data);

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
R<- matrix(c(0,1,0,0,0),nrow=5);
alpha<- 0.05;

#    taus <- c(2:18)/20; 
#    taus <- c(4:36)/40; 
taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  summary(qrfit);
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  #        sigma0tau <- sigma0(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  #        V0 <- solve(jacobtau) %*% sigma0tau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20; 

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);
res0.to.plot<- table0.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);



detach(data);
rm(data);

##########Saving variables for joint graph ################################################################################################################

res.to.plot.80.educ <- res.to.plot;
Kalpha80.educ <- Kalpha;
res0.to.plot.80.educ <- res0.to.plot;
olscoeffs80.educ <- olscoeffs; 


##########################################################################################################################

# Part II:  Implementation (Census 90)

library(foreign);
library(quantreg);
memory.limit(size=500000000);

set.seed(8);

data<-read.dta("census90.dta");
attach(data);

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
R<- matrix(c(0,1,0,0,0),nrow=5);
alpha<- 0.05;

#    taus <- c(2:18)/20; 
#    taus <- c(4:36)/40; 
taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  summary(qrfit);
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  #        sigma0tau <- sigma0(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  #        V0 <- solve(jacobtau) %*% sigma0tau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20; 

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);
res0.to.plot<- table0.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);



detach(data);
rm(data);


##########Saving variables for joint graph ################################################################################################################

res.to.plot.90.educ <- res.to.plot;
Kalpha90.educ <- Kalpha;
res0.to.plot.90.educ <- res0.to.plot;
olscoeffs90.educ <- olscoeffs; 



##########################################################################################################################

# Part II:  Implementation (Census 00)

library(foreign);
library(quantreg);
memory.limit(size=500000000);

set.seed(8);

data<-read.dta("census00.dta");
attach(data);

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
R<- matrix(c(0,1,0,0,0),nrow=5);
alpha<- 0.05;

#    taus <- c(2:18)/20; 
#    taus <- c(4:36)/40; 
taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  summary(qrfit);
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  #        sigma0tau <- sigma0(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  #        V0 <- solve(jacobtau) %*% sigma0tau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20; 

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);
res0.to.plot<- table0.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);



detach(data);
rm(data);


##########Saving variables for joint graph ################################################################################################################

res.to.plot.00.educ <- res.to.plot;
Kalpha00.educ <- Kalpha;
res0.to.plot.00.educ <- res0.to.plot;
olscoeffs00.educ <- olscoeffs; 

##########################################################################################################################

# Part II:  Implementation (Census 10)

library(foreign);
library(quantreg);
memory.limit(size=500000000);

set.seed(8);

data<-read.csv("2010male.csv");
attach(data);
data<- within(data, rm("X"))

n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
R<- matrix(c(0,1,0,0,0),nrow=5);
alpha<- 0.05;

#    taus <- c(2:18)/20; 
#    taus <- c(4:36)/40; 
taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  summary(qrfit);
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  #        sigma0tau <- sigma0(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  #        V0 <- solve(jacobtau) %*% sigma0tau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20; 

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);
res0.to.plot<- table0.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);



detach(data);
rm(data);


##########Saving variables for joint graph ################################################################################################################

res.to.plot.10.educ <- res.to.plot;
Kalpha10.educ <- Kalpha;
res0.to.plot.10.educ <- res0.to.plot;
olscoeffs10.educ <- olscoeffs; 

##########################################################################################################################

# Part II:  Implementation (Census 19)

library(foreign);
library(quantreg);
memory.limit(size=500000000);

set.seed(8);

data<-read.csv("2019male.csv");
attach(data);
data<- within(data, rm("X"))


n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
R<- matrix(c(0,1,0,0,0),nrow=5);
alpha<- 0.05;

#    taus <- c(2:18)/20; 
#    taus <- c(4:36)/40; 
taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  summary(qrfit);
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  #        sigma0tau <- sigma0(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  #        V0 <- solve(jacobtau) %*% sigma0tau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20; 

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);
res0.to.plot<- table0.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);



detach(data);
rm(data);


##########Saving variables for joint graph ################################################################################################################

res.to.plot.19.educ <- res.to.plot;
Kalpha19.educ <- Kalpha;
res0.to.plot.19.educ <- res0.to.plot;
olscoeffs19.educ <- olscoeffs; 




# Part II:  Implementation WOMEN (Census 10)

library(foreign);
library(quantreg);
memory.limit(size=500000000);

set.seed(8);

data<-read.csv("2010female.csv");
attach(data);
data<- within(data, rm("X"))


n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
R<- matrix(c(0,1,0,0,0),nrow=5);
alpha<- 0.05;

#    taus <- c(2:18)/20; 
#    taus <- c(4:36)/40; 
taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  summary(qrfit);
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  #        sigma0tau <- sigma0(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  #        V0 <- solve(jacobtau) %*% sigma0tau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20; 

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);
res0.to.plot<- table0.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);



detach(data);
rm(data);


##########Saving variables for joint graph ################################################################################################################

res.to.plot.fe10.educ <- res.to.plot;
Kalphafe10.educ <- Kalpha;
res0.to.plot.fe10.educ <- res0.to.plot;
olscoeffsfe10.educ <- olscoeffs; 



# Part II:  Implementation WOMEN (Census 19)

library(foreign);
library(quantreg);
memory.limit(size=500000000);

set.seed(8);

data<-read.csv("2019female.csv");
attach(data);
data<- within(data, rm("X"))


n<- dim(data)[1];
B<- 500;
b<- round(5*n^(2/5));
R<- matrix(c(0,1,0,0,0),nrow=5);
alpha<- 0.05;

#    taus <- c(2:18)/20; 
#    taus <- c(4:36)/40; 
taus <- c(10:90)/100;
ntaus <- length(taus);

formula <- logwk~educ+exper+exper2+black;

K <- NULL;

for (i in 1:ntaus) {;
  qrfit <- rq(formula, tau = taus[i], data=data, weights=perwt, method = "fn");
  summary(qrfit);
  coeffs <- coefficients(qrfit);
  res <- residuals(qrfit);
  sigmatau <- sigma(data, n, taus[i], res);
  #        sigma0tau <- sigma0(data, n, taus[i], res);
  jacobtau <- jacobian(data, n, taus[i], res, alpha);
  V <- solve(jacobtau) %*% sigmatau %*% solve(jacobtau);
  #        V0 <- solve(jacobtau) %*% sigma0tau %*% solve(jacobtau);
  K <- rbind(K, subsamplek(data=data, n=n, b=b, B=B, formula, R=R, V, taus[i], coeffs));
};
Kmax <- apply(K,2,FUN=max, na.rm = T);
Kalpha <- quantile(Kmax, 1-alpha);

#ols results

z <- summary(lm(formula, data = data, weights=perwt));
olscoeffs <- t(R) %*% coefficients(z);
vars <- names(z$coef[, 1]);
p <- length(z$coefficients[, 1]);

taus <- c(2:18)/20; 

res.to.plot<- table.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);
res0.to.plot<- table0.rq.res (formula, taus=taus, method="fn", data = data, R=R, alpha, n, hs=F);



detach(data);
rm(data);


##########Saving variables for joint graph ################################################################################################################

res.to.plot.fe19.educ <- res.to.plot;
Kalphafe19.educ <- Kalpha;
res0.to.plot.fe19.educ <- res0.to.plot;
olscoeffsfe19.educ <- olscoeffs; 





############## graph for the intercept and the schooling coefficient ########################################################################################################



postscript("figure2aNEW.eps",horizontal=FALSE,onefile=FALSE,height=8,width=6,pointsize=10);


b80<-  100*res.to.plot.80.educ$tab;
ub80.p<- b80 + 100*Kalpha80.educ*res.to.plot.80.educ$setab;
ub80.m<- b80 - 100*Kalpha80.educ*res.to.plot.80.educ$setab;
b90<-  100*res.to.plot.90.educ$tab;
ub90.p<- b90 + 100*Kalpha90.educ*res.to.plot.90.educ$setab;
ub90.m<- b90 - 100*Kalpha90.educ*res.to.plot.90.educ$setab;
b00<-  100*res.to.plot.00.educ$tab;
ub00.p<- b00 + 100*Kalpha00.educ*res.to.plot.00.educ$setab;
ub00.m<- b00 - 100*Kalpha00.educ*res.to.plot.00.educ$setab;
b10<-  100*res.to.plot.10.educ$tab;
ub10.p<- b10 + 100*Kalpha10.educ*res.to.plot.10.educ$setab;
ub10.m<- b10 - 100*Kalpha10.educ*res.to.plot.10.educ$setab;
b19<-  100*res.to.plot.19.educ$tab;
ub19.p<- b19 + 100*Kalpha19.educ*res.to.plot.19.educ$setab;
ub19.m<- b19 - 100*Kalpha19.educ*res.to.plot.19.educ$setab;
plot( c(0.1,0.9),range(c(6,18)), xlim =c(0.1, 0.9), type="n",xlab="Quantile Index", ylab="Schooling Coefficient (%)", main="");
title(main="                                      SCHOOLING COEFFICIENTS",cex.main=1);
polygon(c(taus,rev(taus)),c(ub80.p,rev(ub80.m)),density=-100, border=F, col='grey60', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub90.p,rev(ub90.m)),density=-100, border=F, col='grey40', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub00.p,rev(ub00.m)),density=-100, border=F, col='grey80', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub10.p,rev(ub10.m)),density=-100, border=F, col='royalblue', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub19.p,rev(ub19.m)),density=-100, border=F, col='royalblue4', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub80.p,rev(ub80.m)),density=50, border=F, col='grey60', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub90.p,rev(ub90.m)),density=50, border=F, col='grey40', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub10.p,rev(ub10.m)),density=50, border=F, col='royalblue', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub19.p,rev(ub19.m)),density=50, border=F, col='royalblue4', lty = 1, lwd = 1);
par(col=1);
lines(taus,b80,lty=6);
lines(taus,b90, lty=1);
lines(taus,b00, lty=2);
lines(taus,b10, lty=2);
lines(taus,b19, lty=1);
abline(h=100*olscoeffs80.educ[1,1],col = 'grey60');
abline(h=100*olscoeffs00.educ[1,1],col = 'grey80');
abline(h=100*olscoeffs90.educ[1,1],col = 'grey40');
abline(h=100*olscoeffs10.educ[1,1],col = 'royalblue');
abline(h=100*olscoeffs19.educ[1,1],col = 'royalblue4');

legend(0.1,18,c("1980","1990","2000","2010","2019"),lwd=c(4,4,4),bty="n",col=c('grey60','grey40','grey80',"royalblue","royalblue4"));
legend(0.1,18,c("1980","1990","2000","2010","2019"),lty=c(6,1,2,2,1),bty="n", lwd = c(1,1,1));




dev.off();

#######women - men comparison
postscript("figure2aWOMEN.eps",horizontal=FALSE,onefile=FALSE,height=8,width=6,pointsize=10);



b10<-  100*res.to.plot.10.educ$tab;
ub10.p<- b10 + 100*Kalpha10.educ*res.to.plot.10.educ$setab;
ub10.m<- b10 - 100*Kalpha10.educ*res.to.plot.10.educ$setab;
b19<-  100*res.to.plot.19.educ$tab;
ub19.p<- b19 + 100*Kalpha19.educ*res.to.plot.19.educ$setab;
ub19.m<- b19 - 100*Kalpha19.educ*res.to.plot.19.educ$setab;
#female
bfe10<-  100*res.to.plot.fe10.educ$tab;
ubfe10.p<- bfe10 + 100*Kalphafe10.educ*res.to.plot.fe10.educ$setab;
ubfe10.m<- bfe10 - 100*Kalphafe10.educ*res.to.plot.fe10.educ$setab;
bfe19<-  100*res.to.plot.fe19.educ$tab;
ubfe19.p<- bfe19 + 100*Kalphafe19.educ*res.to.plot.fe19.educ$setab;
ubfe19.m<- bfe19 - 100*Kalphafe19.educ*res.to.plot.fe19.educ$setab;

plot( c(0.1,0.9),range(c(10,20)), xlim =c(0.1, 0.9), type="n",xlab="Quantile Index", ylab="Schooling Coefficient (%)", main="");
title(main="                                      SCHOOLING COEFFICIENTS",cex.main=1);

polygon(c(taus,rev(taus)),c(ub10.p,rev(ub10.m)),density=-100, border=F, col='royalblue', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub19.p,rev(ub19.m)),density=-100, border=F, col='royalblue4', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ubfe10.p,rev(ubfe10.m)),density=-100, border=F, col='pink', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ubfe19.p,rev(ubfe19.m)),density=-100, border=F, col='pink4', lty = 1, lwd = 1);

polygon(c(taus,rev(taus)),c(ub10.p,rev(ub10.m)),density=50, border=F, col='royalblue', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub19.p,rev(ub19.m)),density=50, border=F, col='royalblue4', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ubfe10.p,rev(ubfe10.m)),density=50, border=F, col='pink', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ubfe19.p,rev(ubfe19.m)),density=50, border=F, col='pink4', lty = 1, lwd = 1);
par(col=1);
lines(taus,b10, lty=2);
lines(taus,b19, lty=2);
lines(taus,bfe10, lty=1);
lines(taus,bfe19, lty=1);

abline(h=100*olscoeffs10.educ[1,1],col = 'royalblue');
abline(h=100*olscoeffs19.educ[1,1],col = 'royalblue4');
abline(h=100*olscoeffsfe10.educ[1,1],col = 'pink');
abline(h=100*olscoeffsfe19.educ[1,1],col = 'pink4');

legend(0.1,20,c("Men 2010","Men 2019","Women 2010", "Women 2019"),lwd=c(4,4,4),bty="n",col=c('royalblue','royalblue4',"pink","pink4"));
legend(0.1,20,c("Men 2010","Men 2019","Women 2010", "Women 2019"),lty=c(2,2,1,1),bty="n", lwd = c(1,1,1));




dev.off();






postscript("figure2bNEW.eps",horizontal=FALSE,onefile=FALSE,height=8,width=6,pointsize=10);



b80<-  res.to.plot.80$tab - res.to.plot.80$tab[10,1];
ub80.p<- b80 + Kalpha.80*res.to.plot.80$setab;
ub80.m<- b80 - Kalpha.80*res.to.plot.80$setab;
b90<-  res.to.plot.90$tab - res.to.plot.90$tab[10,1];
ub90.p<- b90 + Kalpha.90*res.to.plot.90$setab;
ub90.m<- b90 - Kalpha.90*res.to.plot.90$setab;
b00<-  res.to.plot.00$tab - res.to.plot.00$tab[10,1];
ub00.p<- b00 + Kalpha.00*res.to.plot.00$setab;
ub00.m<- b00 - Kalpha.00*res.to.plot.00$setab;
b10<-  res.to.plot.10$tab - res.to.plot.10$tab[10,1];
ub10.p<- b10 + Kalpha.10*res.to.plot.10$setab;
ub10.m<- b10 - Kalpha.10*res.to.plot.10$setab;
b19<-  res.to.plot.19$tab - res.to.plot.19$tab[10,1];
ub19.p<- b19 + Kalpha.19*res.to.plot.19$setab;
ub19.m<- b19 - Kalpha.19*res.to.plot.19$setab;
plot( c(0.1,0.9),range(c(-.8,.7)), xlim =c(0.1, 0.9), yaxp = c(-.8,.8,4), type="n",xlab="Quantile Index", ylab="Log earnings", main="");
title(main="                                                                   CONDITIONAL QUANTILES (at covariate means)",cex.main=1);
polygon(c(taus,rev(taus)),c(ub80.p,rev(ub80.m)),density=-100, border=F, col='grey60', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub90.p,rev(ub90.m)),density=-100, border=F, col='grey40', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub00.p,rev(ub00.m)),density=-100, border=F, col='grey80', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub10.p,rev(ub10.m)),density=-100, border=F, col='royalblue', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub19.p,rev(ub19.m)),density=-100, border=F, col='royalblue4', lty = 1, lwd = 1);

lines(taus,b80,lty=6);
lines(taus,b90,lty=1);
lines(taus,b00,lty=2);
lines(taus,b10,lty=3);
lines(taus,b19,lty=4);

abline(h=0);
legend(0.1,.7,c("1980","1990","2000","2010","2019"),lwd=c(4,4,4),bty="n",col=c('grey60','grey40','grey80',"royalblue","royalblue4"));
legend(0.1,.7,c("1980","1990","2000","2010","2019"),lty=c(6,1,2),bty="n");


dev.off();


postscript("figure2bWOMENNEW.eps",horizontal=FALSE,onefile=FALSE,height=8,width=6,pointsize=10);




b10<-  res.to.plot.10$tab - res.to.plot.10$tab[10,1];
ub10.p<- b10 + Kalpha.10*res.to.plot.10$setab;
ub10.m<- b10 - Kalpha.10*res.to.plot.10$setab;
b19<-  res.to.plot.19$tab - res.to.plot.19$tab[10,1];
ub19.p<- b19 + Kalpha.19*res.to.plot.19$setab;
ub19.m<- b19 - Kalpha.19*res.to.plot.19$setab;
#female
bfe10<-  res.to.plot.fe10$tab - res.to.plot.fe10$tab[10,1];
ubfe10.p<- bfe10 + Kalpha.fe10*res.to.plot.fe10$setab;
ubfe10.m<- bfe10 - Kalpha.fe10*res.to.plot.fe10$setab;
bfe19<-  res.to.plot.fe19$tab - res.to.plot.fe19$tab[10,1];
ubfe19.p<- bfe19 + Kalpha.fe19*res.to.plot.fe19$setab;
ubfe19.m<- bfe19 - Kalpha.fe19*res.to.plot.fe19$setab;
plot( c(0.1,0.9),range(c(-.8,.7)), xlim =c(0.1, 0.9), yaxp = c(-.8,.8,4), type="n",xlab="Quantile Index", ylab="Log earnings", main="");
title(main="                                                                   CONDITIONAL QUANTILES (at covariate means)",cex.main=1);

polygon(c(taus,rev(taus)),c(ub10.p,rev(ub10.m)),density=-100, border=F, col='royalblue', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ub19.p,rev(ub19.m)),density=-100, border=F, col='royalblue4', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ubfe10.p,rev(ubfe10.m)),density=-100, border=F, col='pink', lty = 1, lwd = 1);
polygon(c(taus,rev(taus)),c(ubfe19.p,rev(ubfe19.m)),density=-100, border=F, col='pink4', lty = 1, lwd = 1);


lines(taus,b10,lty=2);
lines(taus,b19,lty=2);
lines(taus,bfe10,lty=1);
lines(taus,bfe19,lty=1);

abline(h=0);
legend(0.1,.7,c("Men 2010","Men 2019","Women 2010", "Women 2019"),lwd=c(4,4,4),bty="n",col=c('royalblue','royalblue4',"pink","pink4"));
legend(0.1,.7,c("Men 2010","Men 2019","Women 2010", "Women 2019"),lty=c(2,2,1,1),bty="n");


dev.off();


