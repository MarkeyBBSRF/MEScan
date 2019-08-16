####### find the parameter a,b of the prior of q_i by maximizing the marginal likelihood over (a,b) ###################################
#srda = "pan12/mutationtable.pan12.Rdata"
source("MADGiC.R")
default = T
if(default){
	####### find the parameter a,b of the prior of q_i by maximizing the marginal likelihood over (a,b) ###################################

	y= try(optim(c(-25,-20),find.prior.param,muttable=muttable,uniqueA=uniqueA,tableA=tableA,S=S), silent=TRUE)
	if(class(y) == "try-error"){ y= optim(c(-20,-10),find.prior.param,muttable=muttable,uniqueA=uniqueA,tableA=tableA,S=S) }
	a=y$par[1]
	b=y$par[2]


	######################################################################################################################################

	sample.name=unique(c(sil.mutab[,8],nonsil.mutab[,8]))     ### sample id
	if(S>length(sample.name))
	  sample.name=c(sample.name,1:(S-length(sample.name)))

	##### calculate bij  ##########################################################################################################
	bijAndQvec <- calculate.bij( gid, gene.rep.expr, p, s, epsilon, delta, exome.constants.mult, sample.name,muttable,nonsil.mut.type,a,b,uniqueA,tableA)
	names(bijAndQvec$bij) <- gid

	save(bijAndQvec,sample.name,gid,file=file.path(path,"Bij_default.Rdata"))

	}else{
		mmin=-3
	y= try(optim(c(-25,-20),find.prior.param,muttable=muttable,uniqueA=uniqueA,tableA=tableA,S=S, M=mmin), silent=TRUE)
	# if(class(y) == "try-error"){ y= optim(c(-20,-10),find.prior.param,muttable=muttable,uniqueA=uniqueA,tableA=tableA,S=S) }
	if(!exists("y")|| class(y) == "try-error" ){
	ttry = c(-30,-10)
	mmin=-5
	while(!exists("y") || class(y)=="try-error" ) {
	   # ttry[1] = ttry[1]-5
	    #print(ttry)
	    #mmin = mmin-2
	    mmin=mmin -1
	    stopifnot(mmin > min(ttry))
	    print(mmin)
	    y= try(optim(ttry,find.prior.param,muttable=muttable,uniqueA=uniqueA,tableA=tableA,S=S,M=mmin), silent=TRUE)
	}
	}
	a=y$par[1]
	b=y$par[2]
	save(y,file = file.path(path,'y.rda'))
	######################################################################################################################################

	sample.name=unique(c(sil.mutab[,8],nonsil.mutab[,8]))     ### sample id
	if(S>length(sample.name))
	  sample.name=c(sample.name,1:(S-length(sample.name)))

	##### calculate bij  ##########################################################################################################
	bijAndQvec <- calculate.bij( gid, gene.rep.expr, p, s, epsilon, delta, exome.constants.mult, sample.name,muttable,nonsil.mut.type,a,b,uniqueA,tableA)
	## check problematic patient
	##unique(unlist(sapply(gene.idindex,function(x) if(length(x)!=5554) return(x))))

	names(bijAndQvec$bij) <- gid
	save(bijAndQvec,sample.name,gid,file=file.path(path,"Bij.Rdata"))


}
