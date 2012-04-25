
SKAT_TruncatedNormal_Get_Moment = function(c1, c2, sigma, Xa){
  a1<-(c1-Xa)/sigma
  a2<-(c2-Xa)/sigma	
	f1<-dnorm(a1)
  f2<-dnorm(a2)
	F1<-pnorm(a1)
	F2<-pnorm(a2)
	f12<-f2-f1
	F12<-F2 + 1 - F1
	M<- f12/F12 * sigma 
	V<- (a2*f2-a1*f1)/F12 + M^2 /sigma^2
	W<- (1 - V) * sigma^2
	Mu<- Xa - M
	return(list(M=M, V=V, W=W, Mu=Mu))
}

#
#	The latest version of the TruncatedNormalSove function
#
SKAT_TruncatedNormalSolve = function(formula, c1, c2, delta=0.001, gamma0="default", MAXITERNUM=20, data=NULL){
	preerror = 0
	X<-model.matrix(formula,data=data)
	Y<-model.frame(formula, data=data)[,1]
	if(is.character(gamma0) && gamma0=="default"){
		out.lm<-lm(formula, data=data)
    sigma.ols =mean(abs(out.lm$residuals))
    c = (c1-c2)/2 
		gamma0 = c(out.lm$coefficients,sigma.ols)
	}
	n = nrow(X)
	K = ncol(X)
	iternum = 0
  sigma = gamma0[K+1]
  alpha0 = gamma0[-(K+1)]
  while(T){
  		Xa = (X%*%alpha0)[,1]
  		J = matrix(0,nrow=K,ncol=K)
  		a1<-(c1-Xa)/sigma
  		a2<-(c2-Xa)/sigma
  		f1<-dnorm(a1)
  		f2<-dnorm(a2)
  		F1<-pnorm(a1)
  		F2<-pnorm(a2)
  		Y_Xa<-(Y-Xa)
  		f12<-f1-f2
  		F12<-F2 + 1 - F1
  		J[1:K, 1:K]<-t(X) %*% ( X * (-1+(a2*f2-a1*f1)/F12+(f12/F12)^2)) / sigma^2
  		Jinv = solve(J)
  		V = matrix(0,nrow=K,ncol=1)
  		V[1:K,1] = colSums((X/sigma)*(Y_Xa/sigma-f12/F12))
  		sigma_next = uniroot(function(sigmaT) sum(-1/sigmaT+(Y_Xa)^2/sigmaT^3+((a2/sigmaT)*f2-((c1-X%*%alpha0)/sigmaT^2)*f1)/F12),c(sigma/20,10*sigma))$root
      alphaNext = alpha0-Jinv%*%V
  		iternum = iternum + 1
  		curerror = sum(abs(alphaNext-alpha0))+abs(sigma-sigma_next)
 		  if(iternum > MAXITERNUM && curerror > preerror){stop("Newton-Raphson diverging. Try different initial guess for gamma0.")}
  		preerror = curerror  
  		if(curerror < delta){
          end = TRUE
          if(end){
	    		  MLEs = c(alphaNext,sigma)
	    		  MLEs = data.frame(MLEs)
	    		  rownames = rep(0,K+1)
	    		  for(i in 1:K){
	      		  rownames[i] = paste("alpha",i-1,sep="")
	    		  }
	    		  rownames[1] = "intercept"
	    		  rownames[K+1] = "sigma"
	    		  row.names(MLEs) = rownames
  			    alpha0 = alphaNext
  			    Xa = (X%*%alpha0)[,1]
			
			  re<-SKAT_TruncatedNormal_Get_Moment(c1, c2, sigma, Xa)
	    		  return(list(coef=MLEs, M = re$M, V= re$V, W=re$W, Mu=re$Mu, Y=Y, X1=X, sigma=sigma, Xa=Xa))
          }
  		}
      alpha0 = alphaNext
      sigma=sigma_next
	}
}

#	Null model (ECP)
#

SKAT_Null_Model_EP = function(formula, c1, c2, delta=0.001, gamma0="default", MAXITERNUM=20 , data=NULL, n.Resampling=0){
	
	
	out<-SKAT_TruncatedNormalSolve(formula, c1, c2, delta, gamma0, MAXITERNUM,data)
	
	res = out$Y - out$Mu
	pi_1 = out$W
	X1 = out$X1
	n<-length(res)
	id_include=1:n

	res.out=NULL
	type.Resampling=NULL
	

	if(n.Resampling > 0){
		res.out<-SKAT_Get_TruncatedNormal(length(res), out$Xa, out$sigma, c1, c2, nSet=n.Resampling) - out$Mu

		res.out<-t(t(res.out) - colMeans(res.out))
	}

  	re<-list(res=res, X1=X1,res.out=res.out,out_type="D", 
	n.Resampling=n.Resampling, type.Resampling=NULL, id_include=id_include, mu=out$Mu,pi_1=pi_1, coef = out$coef, Xa=out$Xa)

	class(re)<-"SKAT_NULL_Model"

	return(re)
	
}


#	Null model with small sample adjustment
#

SKAT_Null_Model_EP_MomentAdjust = function(formula, c1, c2, delta=0.001, gamma0="default", MAXITERNUM=20 , data=NULL, n.Resampling=0){
	
	
	out<-SKAT_TruncatedNormalSolve(formula, c1, c2, delta, gamma0, MAXITERNUM,data)
	
	res = out$Y - out$Mu
	pi_1 = out$W
	X1 = out$X1
	n<-length(res)
	id_include=1:n

	res.out=NULL
	type.Resampling=NULL
	

	if(n.Resampling > 0){
		res.out<-SKAT_Get_TruncatedNormal(length(res), out$Xa, out$sigma, c1, c2, nSet=n.Resampling) - out$Mu

		res.out<-t(t(res.out) - colMeans(res.out))
	}

  	re1<-list(res=res, X1=X1,res.out=res.out,out_type="D", 
	n.Resampling=n.Resampling, type.Resampling=NULL, id_include=id_include, mu=out$Mu,pi_1=pi_1, coef = out$coef, Xa=out$Xa)


	res.out1<-SKAT_Get_TruncatedNormal(length(res), out$Xa, out$sigma, c1, c2, nSet=5000) - out$Mu
	res.out1<-t(t(res.out1) - colMeans(res.out1))

	re2<-list(res.out=res.out1)
	re<-list(re1=re1, re2=re2, is_kurtosis_adj= TRUE, type = "ECP")


	class(re)<-"SKAT_NULL_Model_ADJ"

	return(re)
	
}


#	Function to simulate truncated normal random variables
#	x > c1 or x < c2

SKAT_Get_TruncatedNormal<-function(n, mu, sigma, c1, c2, nSet, nUpper = n/2){


	prob1<-1- pnorm(c1,mean=mu, sd=sigma) 
	prob2<-pnorm(c2,mean=mu, sd=sigma) 
	
	prob<-prob1/(prob1 + prob2)

	temp.b<-rbinom(n*nSet,1,rep(prob,nSet))	# Prob x > c1

	lower = rep(- 10^5,n*nSet) * (1-temp.b) + c1
	upper = rep(10^5,n*nSet) * temp.b + c2
	
	re<-rtnorm(n*nSet, mean=mu, sd=sigma, lower=lower, upper=upper)
	re.M<-matrix(re, ncol=nSet, byrow=FALSE)

	cutoff<<-cbind(lower,upper, temp.b, mu)
	rm(re)
	rm(lower)
	rm(upper)
	
	return(re.M)

}


## Rejection sampling algorithm by Robert (Stat. Comp (1995), 5, 121-5)
## for simulating from the truncated normal distribution.
## This code come from msm package.

rtnorm <- function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
    if (length(n) > 1)
        n <- length(n)
    mean <- rep(mean, length=n)
    sd <- rep(sd, length=n)
    lower <- rep(lower, length=n)
    upper <- rep(upper, length=n)
    lower <- (lower - mean) / sd ## Algorithm works on mean 0, sd 1 scale
    upper <- (upper - mean) / sd
    ind <- seq(length=n)
    ret <- numeric(n)
    ## Different algorithms depending on where upper/lower limits lie.
    alg <- ifelse(
                  lower > upper,
                  -1,# return NaN if lower > upper
                  ifelse(
                         ((lower < 0 & upper == Inf) |
                          (lower == -Inf & upper > 0) |
                          (is.finite(lower) & is.finite(upper) & (lower < 0) & (upper > 0) & (upper-lower > sqrt(2*pi)))
                          ),
                         0, # standard "simulate from normal and reject if outside limits" method. Use if bounds are wide.  FIXME HSOULD BE 
                         ifelse(
                                (lower >= 0 & (upper > lower + 2*sqrt(exp(1)) /
                                 (lower + sqrt(lower^2 + 4)) * exp((lower*2 - lower*sqrt(lower^2 + 4)) / 4))),
                                1, # rejection sampling with exponential proposal. Use if lower >> mean
                                ifelse(upper <= 0 & (-lower > -upper + 2*sqrt(exp(1)) /
                                       (-upper + sqrt(upper^2 + 4)) * exp((upper*2 - -upper*sqrt(upper^2 + 4)) / 4)),
                                       2, # rejection sampling with exponential proposal. Use if upper << mean.
                                       3)))) # rejection sampling with uniform proposal. Use if bounds are narrow and central.
    
    ind.nan <- ind[alg==-1]; ind.no <- ind[alg==0]; ind.expl <- ind[alg==1]; ind.expu <- ind[alg==2]; ind.u <- ind[alg==3]
    ret[ind.nan] <- NaN
    while (length(ind.no) > 0) {
        y <- rnorm(length(ind.no))
        done <- which(y >= lower[ind.no] & y <= upper[ind.no])
        ret[ind.no[done]] <- y[done]
        ind.no <- setdiff(ind.no, ind.no[done])
    }
    stopifnot(length(ind.no) == 0)
    while (length(ind.expl) > 0) {
        a <- (lower[ind.expl] + sqrt(lower[ind.expl]^2 + 4)) / 2
        z <- rexp(length(ind.expl), a) + lower[ind.expl]
        u <- runif(length(ind.expl))
        done <- which((u <= exp(-(z - a)^2 / 2)) & (z <= upper[ind.expl]))
        ret[ind.expl[done]] <- z[done]
        ind.expl <- setdiff(ind.expl, ind.expl[done])
    }
    stopifnot(length(ind.expl) == 0)
    while (length(ind.expu) > 0) {
        a <- (-upper[ind.expu] + sqrt(upper[ind.expu]^2 +4)) / 2
        z <- rexp(length(ind.expu), a) - upper[ind.expu]
        u <- runif(length(ind.expu))
        done <- which((u <= exp(-(z - a)^2 / 2)) & (z <= -lower[ind.expu]))
        ret[ind.expu[done]] <- -z[done]
        ind.expu <- setdiff(ind.expu, ind.expu[done])
    }
    stopifnot(length(ind.expu) == 0)
    while (length(ind.u) > 0) {
        z <- runif(length(ind.u), lower[ind.u], upper[ind.u])
        rho <- ifelse(lower[ind.u] > 0,
                      exp((lower[ind.u]^2 - z^2) / 2), ifelse(upper[ind.u] < 0,
                                                            exp((upper[ind.u]^2 - z^2) / 2),
                                                            exp(-z^2/2)))
        u <- runif(length(ind.u))
        done <- which(u <= rho)
        ret[ind.u[done]] <- z[done]
        ind.u <- setdiff(ind.u, ind.u[done])
    }
    stopifnot(length(ind.u) == 0)
    ret*sd + mean
}


