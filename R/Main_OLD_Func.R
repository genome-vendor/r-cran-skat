SKAT_OLD = function(Z,y, X=NULL, kernel = "linear.weighted", out_type="C",method="davies", weights.beta=c(1,25) , weights = NULL, is_intercept = TRUE){

	re<-SKAT_MAIN(Z,y, X, kernel, out_type,method, weights.beta, weights,is_intercept=is_intercept)	
	return(re)
}


#
#	4 methods (keep in mind)
#
SKAT_MAIN = function(Z,y, X=NULL, kernel = "linear.weighted", out_type="C",method="davies", weights.beta=c(1,25) , weights = NULL, impute.method = "fixed", SetID = NULL, is_intercept = TRUE, r.corr=0, is_check_genotype = TRUE, is_dosage = FALSE){

	warning("It is an old interface!, please run SKAT with SKAT_Null_Model or SKAT_Null_Model_MomentAdjust!")

	out.method<-SKAT_Check_Method(method,r.corr)
	method=out.method$method
	r.corr=out.method$r.corr

	SKAT_Check_RCorr(kernel, r.corr)
	SKAT_MAIN_Check_OutType(out_type)
	out.xy<-SKAT_MAIN_Check_XY(y,X,is_intercept)
	out.z<-SKAT_MAIN_Check_Z(Z,out.xy$n,out.xy$id_missxy
	,out.xy$id_include, SetID,weights,weights.beta, impute.method, is_check_genotype, is_dosage)

	if(out.z$return ==1){
		return(out.z)
	}

	if(out_type == "C"){
		re<-SKAT.linear(out.z$Z.test,out.xy$y.test,out.xy$X1.test,kernel = kernel, weights = out.z$weights, method=method, r.corr=r.corr)
	} else if (out_type == "D"){
		re<-SKAT.logistic(out.z$Z.test,out.xy$y.test,out.xy$X1.test, kernel = kernel, weights = out.z$weights, method=method, r.corr=r.corr)
	}

	return(re)

}


#
#	Check the parameter y and X
#
#
SKAT_MAIN_Check_XY<-function(y,X,is_intercept){

	n = length(y)
	if (!is.null(X)) {
    		if (class(X)!= "matrix") stop("X is not a matrix")
		if(sum((X[,1] - 1)^2,na.rm = TRUE) == 0){
			X1 = X
		} else if(is_intercept){
    			X1 = cbind(1,X)
		} else {
    			X1 = X
		}

  	} else if(is_intercept == TRUE){
    		X1 = as.matrix(rep(1,n))
		
  	} else {
    		stop("Error: Additional covariates (X) are needed when is_intercept = FALSE")
  	}

	if (nrow(X1)!=n) stop("Dimensions of y and X do not match")

	# Check missing of y and X
	id_missy<-which(is.na(y))
	id_missx<-NULL
	for(k in 1:dim(X1)[2]){
		temp1<-which(is.na(X1[,k]))
		id_missx<-union(id_missx,temp1)
	}
	id_missxy<-sort(union(id_missx,id_missy))
	
	id_include<-1:n
	if(length(id_missxy) > 0){
		MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",length(id_missxy))
		warning(MSG,call.=FALSE)
		id_include<-id_include[-id_missxy]
		
	}

	y.test<-y[id_include]
	X1.test<-as.matrix(X1[id_include,])

	return(list(n=n,X1=X1,id_missxy=id_missxy,id_include=id_include,
	y.test=y.test,X1.test=X1.test))

}

