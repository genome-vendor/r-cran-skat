
SKAT_Null_Model_SmallSample = function(formula, data=NULL){
	
	
	# check missing 
	obj1<-model.frame(formula,na.action = na.omit,data=NULL)
	obj2<-model.frame(formula,na.action = na.pass,data=NULL)

	n<-dim(obj2)[1]
	n1<-dim(obj1)[1]
	id_include<-as.numeric(rownames(obj1))

	if(n - n1 > 0){
		MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n - n1)
		warning(MSG,call.=FALSE)
	}

	is_skewness_adj == TRUE
	if(n > 2000){
		is_skewness_adj == FALSE
	} 

	re1<-Get_SKAT_Residuals.logistic (formula, data, n.Resampling, type.Resampling, id_include )
	re2<-NULL

	if(is_skewness_adj == TRUE){
		re2<-Get_SKAT_Residuals.logistic (formula, data, 10000, type.Resampling, id_include )
	}

	re<-list(re1=re1, re2=re2, is_skewness_adj= is_skewness_adj)

	class(re)<-"SKAT_NULL_Model_SMALL"
	return(re)
	
}


SKAT_With_NullModel_SMALL = function(Z, obj.res.a, kernel = "linear.weighted", method="davies", weights.beta=c(1,25), weights = NULL,
impute.method = "random", SetID = NULL, r.corr=0, is_check_genotype=TRUE){

	re<-SKAT_With_NullModel_ADJ(Z, obj.res.a, kernel, method, weights.beta, weights,
impute.method , SetID , r.corr, is_check_genotype)

	#out<-list(p.value=p.value,  p.value.noadj= p.value.noadj
	#, p.value.d=out.d1$p.value, p.value.d.noadj=out.d2$p.value
	#, muQ=muQ, varQ=varQ, df=df)
	
}
 


