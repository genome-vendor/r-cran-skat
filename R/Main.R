
SKAT = function(Z,obj,...){

	ml<-match.call()
	ml.name<-names(ml)
	IDX1<-which(ml.name == "y")
	if(length(IDX1) > 0){
		re<-SKAT_MAIN(Z,...)
	} else {

		if(class(obj) == "SKAT_NULL_Model_ADJ"){
			re<-SKAT_With_NullModel_ADJ(Z,obj, ...)
		} else if(class(obj) == "SKAT_NULL_Model"){
			re<-SKAT_With_NullModel(Z,obj, ...)
		} else {
			#re<-SKAT_MAIN(Z,obj, ...)
			stop("The old interface is defunct! Please run SKAT_NULL_Model first!")
		}

	}	
	class(re)<-"SKAT_OUT"
	return(re)
}

#
#	Check the out_type
#
SKAT_MAIN_Check_OutType<-function(out_type){
 	
	if(out_type != "C" && out_type != "D"){
		stop("Invalid out_type!. Please use either \"C\" for the continous outcome or \"D\" for the dichotomous outcome.")
	}

}

#
#	Check the Z, and do imputation
#
#
SKAT_MAIN_Check_Z<-function(Z, n, id_include, SetID, weights, weights.beta, impute.method, is_check_genotype, is_dosage, missing_cutoff){

	#############################################
	# Check parameters

	if (class(Z)!= "matrix") stop("Z is not a matrix")
	if (nrow(Z)!=n) stop("Dimensions of y and Z do not match")
 	if(is_dosage ==TRUE){
		impute.method="fixed"
	}
	#####################################################
	# Check Z

	if(!is_check_genotype && !is_dosage){
		Z.test<-Z[id_include,]
		if(!is.matrix(Z.test)){
			Z.test<-as.matrix(Z.test)
		}
		return(list(Z.test=Z.test,weights=weights, return=0) )
	}

	##############################################
	# Check Missing 

	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 

	###################################################
	# Check missing rates and exclude any SNPs with missing rate > missing_cutoff
	# Also exclude non-polymorphic SNPs
	m = ncol(Z)
	ID_INCLUDE_SNP<-NULL
	for(i in 1:m){
		missing.ratio<-length(which(is.na(Z[,i])))/n
		sd1<-sd(Z[,i], na.rm=TRUE)
		if(missing.ratio < missing_cutoff && sd1 > 0){
			ID_INCLUDE_SNP<-c(ID_INCLUDE_SNP,i)
		}
	}
	
	if(length(ID_INCLUDE_SNP) == 0){

		if(is.null(SetID)){
			msg<-sprintf("ALL SNPs have either high missing rates or no-variation. P-value=1")
		} else {
			msg<-sprintf("In %s, ALL SNPs have either high missing rates or no-variation. P-value=1",SetID )
		}
		warning(msg,call.=FALSE)
		
		re<-list(p.value = 1, p.value.resampling =NA, Test.Type = NA, Q = NA, param=list(n.marker=0, n.marker.test=0), return=1 )   

	} else if(m - length(ID_INCLUDE_SNP) > 0 ){

		if(is.null(SetID)){
			msg<-sprintf("%d SNPs with either high missing rates or no-variation are excluded!", m - length(ID_INCLUDE_SNP))
		} else {
			msg<-sprintf("In %s, %d SNPs with either high missing rates or no-variation are excluded!",SetID, m - length(ID_INCLUDE_SNP) )
		}

		warning(msg,call.=FALSE)	
		Z<-as.matrix(Z[,ID_INCLUDE_SNP])
	}


	##################################################################
	# doing imputation

	MAF<-colMeans(Z, na.rm = TRUE)/2
	MAF1<-colMeans(as.matrix(Z[id_include,]),na.rm=TRUE)/2
	IDX.Err<-which(MAF > 0.5)	
	if(length(IDX.Err) > 0){
		#msg<-sprintf("Genotypes of some variants are not the number of minor allele! It is fixed!")
		msg<-sprintf("Genotypes of some variants are not the number of minor alleles!")
		warning(msg,call.=FALSE)

		# Fixed by SLEE
		#Z[,IDX.Err]<-2 - Z[,IDX.Err]
		#MAF[IDX.Err]<-1- MAF[IDX.Err]
	}

	###########################################
	# Check non-polymorphic

	if(length(which(MAF1 > 0)) == 0){
		
		if(is.null(SetID)){
			msg<-sprintf("No polymorphic SNP. P-value = 1" )
		} else {
			msg<-sprintf("In %s, No polymorphic SNP. P-value = 1",SetID )
		}
		warning(msg,call.=FALSE)
		re<-list(p.value = 1, p.value.resampling =NA, Test.Type = NA, Q = NA, param=list(n.marker=0, n.marker.test=0), return=1 )   
		return(re)
	}

	##########################################
	# Missing Imputation
	if(length(IDX_MISS) > 0){
		if(is.null(SetID)){
			msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )
		} else {
			msg<-sprintf("In %s, the missing genotype rate is %f. Imputation is applied.", SetID, (length(IDX_MISS))/length(Z) )
		}

		warning(msg,call.=FALSE)
		Z<-Impute(Z,impute.method)
	} 
	
	##########################################
	# Get Weights

	if(is.null(weights)){
		weights<-Beta.Weights(MAF,weights.beta)
	}

	###########################################
	# Check missing of y and X

	if(n - length(id_include)  > 0){
	
		id_Z<-which(MAF1 > 0)

		if(length(id_Z) == 0){

			if(is.null(SetID)){
				msg<-sprintf("No polymorphic SNP. P-value = 1" )
			} else {
				msg<-sprintf("In %s, No polymorphic SNP. P-value = 1",SetID )
			}
			warning(msg,call.=FALSE)
			re<-list(p.value = 1, p.value.resampling =NA, Test.Type = NA, Q = NA, param=list(n.marker=0, n.marker.test=0), return=1 )   

		} else if (length(id_Z) == 1){
			Z<-cbind(Z[,id_Z])
		} else {
			Z<-Z[,id_Z]
		}

		if(!is.null(weights)){
			weights<-weights[id_Z]
		}

	}	
	
	if( dim(Z)[2] == 1){

		if(is.null(SetID)){
			msg<-sprintf("Only one SNP in the SNP set!" )
		} else {
			msg<-sprintf("In %s, Only one SNP in the SNP set!"
			,SetID )
		}
		warning(msg,call.=FALSE)

		Z.test<-as.matrix(Z[id_include,])

	} else {

		Z.test<-Z[id_include,]

	}

	return(list(Z.test=Z.test,weights=weights, return=0) )

}

SKAT_Check_RCorr<-function(kernel, r.corr){

	if(length(r.corr) == 1 && r.corr[1] == 0){
		return(1)
	}
	if(kernel != "linear" && kernel != "linear.weighted"){
		stop("Error: non-zero r.corr only can be used with linear or linear.weighted kernels")
	}

	for(i in 1:length(r.corr)){
		if(r.corr[i] < 0 || r.corr[i] > 1){
			stop("Error: r.corr should be >= 0 and <= 1")
		}
	}



}

SKAT_Check_Method<-function(method,r.corr){


	if(method != "liu"  && method != "davies" && method != "liu.mod" && method != "optimal" && method != "optimal.moment" && method != "adjust"  ){
		stop("Invalid method!")
	}
	
	if((method == "optimal" || method =="optimal.moment") && length(r.corr) == 1){
		r.corr = (0:10)/10
		#r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)
	}	
	if(method =="optimal"){
		method="davies"
	} else if (method =="optimal.moment") {
		method="liu.mod"
	}

	re<-list(method=method,r.corr=r.corr)
	return(re)

}




SKAT_With_NullModel = function(Z, obj.res, kernel = "linear.weighted", method="davies", weights.beta=c(1,25), weights = NULL,
impute.method = "fixed", SetID = NULL, r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15){

	
	n<-dim(Z)[1]
	m<-dim(Z)[2]

	out.method<-SKAT_Check_Method(method,r.corr)
	method=out.method$method
	r.corr=out.method$r.corr


	SKAT_Check_RCorr(kernel, r.corr)

	out.z<-SKAT_MAIN_Check_Z(Z, n, obj.res$id_include, SetID, weights, weights.beta, impute.method, is_check_genotype, is_dosage, missing_cutoff)
	if(out.z$return ==1){
		out.z$param$n.marker<-m
		return(out.z)
	}

	if(length(r.corr) > 1 && dim(out.z$Z.test)[2] <= 1){
		r.corr=0
	}

	if(obj.res$out_type == "C"){
		  if( (kernel =="linear" || kernel == "linear.weighted") && n > m){
		    re = SKAT.linear.Linear(obj.res$res,out.z$Z.test
			,obj.res$X1, kernel, out.z$weights,obj.res$s2,method
			,obj.res$res.out, obj.res$n.Resampling,r.corr=r.corr)
		  } else {  
		    re = SKAT.linear.Other(obj.res$res,out.z$Z.test
			,obj.res$X1, kernel, out.z$weights,obj.res$s2,method
			,obj.res$res.out, obj.res$n.Resampling)  
		  }
	} else if (obj.res$out_type == "D"){

		if( (kernel =="linear" || kernel == "linear.weighted") && n > m){
			re = SKAT.logistic.Linear(obj.res$res, out.z$Z.test
			,obj.res$X1, kernel, out.z$weights, obj.res$pi_1,method
			,obj.res$res.out, obj.res$n.Resampling,r.corr=r.corr)
		} else {  
			re = SKAT.logistic.Other(obj.res$res,out.z$Z.test
			,obj.res$X1, kernel, out.z$weights, obj.res$pi_1, method
			,obj.res$res.out, obj.res$n.Resampling)  
		}
	}

	re$param$n.marker<-m
	re$param$n.marker.test<-dim(out.z$Z.test)[2]
	return(re)

}
 
#
#	Adjustment methods only use liu.mod, so it doesn't need method the "method" parameter
#	I use this field for outcome.type for subfunctions
#	
SKAT_With_NullModel_ADJ = function(Z, obj.res.a, kernel = "linear.weighted", method="adjust", weights.beta=c(1,25), weights = NULL,
impute.method = "fixed", SetID = NULL, r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15){

	
	n<-dim(Z)[1]
	m<-dim(Z)[2]
	obj.res<-obj.res.a$re1

	out.method<-SKAT_Check_Method(method,r.corr)
	method=out.method$method
	r.corr=out.method$r.corr

	SKAT_Check_RCorr(kernel, r.corr)
	# Use method field for the type of outcome
	method = obj.res.a$type

	out.z<-SKAT_MAIN_Check_Z(Z, n, obj.res$id_include, SetID, weights, weights.beta, impute.method, is_check_genotype, is_dosage, missing_cutoff)
	if(out.z$return ==1){
		out.z$param$n.marker<-m
		return(out.z)
	}

	res2<-NULL
	if(obj.res.a$is_kurtosis_adj){
		res2<-obj.res.a$re2$res.out
	}	

	if(length(r.corr) > 1 && dim(out.z$Z.test)[2] <= 1){
		r.corr=0
	}

	if(length(r.corr) == 1){

		re = KMTest.logistic.Linear.VarMatching (obj.res$res,out.z$Z.test
			, obj.res$X1, kernel, out.z$weights, obj.res$pi_1, method
			, obj.res$res.out, obj.res$n.Resampling,r.corr=r.corr
			, obj.res$mu, res.moments=res2)

	} else {

		re = SKAT_Optimal_Logistic_VarMatching(obj.res$res, out.z$Z.test
			, obj.res$X1, kernel, out.z$weights, obj.res$pi_1, method
			, obj.res$res.out, obj.res$n.Resampling, r.corr, obj.res$mu
			, res.moments=res2)

	}

	re$param$n.marker<-m
	re$param$n.marker.test<-dim(out.z$Z.test)[2]
	return(re)


}
 
 
