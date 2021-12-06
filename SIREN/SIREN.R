BMatrixFn=function(expMatrix)
{
cat("calculating B-matrix...\n")
library(splines)
nc=ncol(expMatrix)
ng=nrow(expMatrix)
bMatrix=array(0,dim=c(nc,10,ng))

for(counterGene in seq(ng))
{
	j=bs(scale(expMatrix[counterGene,]),df=10,degree=2)
	for(counterCondition in seq(nc))
	{
		bMatrix[counterCondition,1,counterGene]=j[counterCondition,1]
		bMatrix[counterCondition,2,counterGene]=j[counterCondition,2]
		bMatrix[counterCondition,3,counterGene]=j[counterCondition,3]
		bMatrix[counterCondition,4,counterGene]=j[counterCondition,4]
		bMatrix[counterCondition,5,counterGene]=j[counterCondition,5]
		bMatrix[counterCondition,6,counterGene]=j[counterCondition,6]
		bMatrix[counterCondition,7,counterGene]=j[counterCondition,7]
		bMatrix[counterCondition,8,counterGene]=j[counterCondition,8]
		bMatrix[counterCondition,9,counterGene]=j[counterCondition,9]
		bMatrix[counterCondition,10,counterGene]=j[counterCondition,10]
	}
}

return(bMatrix)
}

PaCalculate=function(bMatrix)
{
cat("Calculating marginal probabilities...\n")
#nc: is the number of conditions
nc=length(bMatrix[,1,1])

prMatrix=matrix(0,nrow=length(bMatrix[1,1,]),ncol=10)
#10: the number of bins for zero and non-zero values
for(rowNumber in seq(length(bMatrix[1,1,])))
for (m in seq(10))   
{

    for(k in seq(nc))   
        prMatrix[rowNumber,m]=prMatrix[rowNumber,m]+bMatrix[k,m,rowNumber]
    
    prMatrix[rowNumber,m]=prMatrix[rowNumber,m]/nc
    
}
return(prMatrix)
}

PabCalculate=function(bMatrixA,bMatrixB)
{

#nc=number of conditions
nc=length(bMatrixA[,1]);

#10: the number of bins for zero and non-zero values
#prAB: joint probability of genes A and B
prAB=matrix(0,ncol=10,nrow=10)

for (m in seq(10))   
{
    for (mm in seq(10))
    {    
        
        for (k in seq(nc))
            prAB[m,mm]=prAB[m,mm]+bMatrixA[k,m]*bMatrixB[k,mm]
        

        prAB[m,mm]=prAB[m,mm]/nc
    }
}

return(prAB)
}

ScoreFn=function(bMatrix,weightMatrix,prMatrix,net=NULL)
{
cat("Calculating SIREN scores...\n")
#Aindex= is the row number of the gene A of interest for calculating it's
#probability

#Bindex= is the row number of the gene B of interest for calculating it's
#probability

#prA: Marginal Probabilities of gene A (10 bins => 10 marginal probabilities)
#prB: Marginal probabilities of gene B (10 bins => 10 marginal probabilities)
#prAB: Joint probabilities of gene A and B (10 bins => 100 (10*10) joint probabilities)

SIRENscore=c(0)

counter=0

if(is.null(net))
{
for(Aindex in seq(length(bMatrix[1,1,])))
{
cat(Aindex)
cat("\n")
	for(Bindex in seq(Aindex+1,length(bMatrix[1,1,]),1))
	{
		counter=counter+1
		SIRENscore[counter]=0
		#SIRENscore[counter,1]=Aindex
		#SIRENscore[counter,2]=Bindex
		#SIRENscore[counter,4]=1
		prAB=PabCalculate(bMatrix[,,Aindex],bMatrix[,,Bindex])
		prA=prMatrix[Aindex,]

		prB=prMatrix[Bindex,]

		#Calculating Score for genes A and B
		#10: the number of bins for zero and non-zero values
		for (m in seq(10))
		{   
		    for (mm in seq(10))
	            {
			xN=log((prAB[m,mm])/(prA[m]*prB[mm]))/(-1*log(prAB[m,mm]))
			tryCatch({
				if(xN>0)
		       		{
					 SIRENscore[counter]=SIRENscore[counter]+(prAB[m,mm])*weightMatrix[m,mm]*xN
				}
			 }
			 ,error=function(e) {}
			 ,finally={})
		    }
		}

	}
}
}
else
{
cat("Network file is provided...")
cat("\n")
for( i in seq(length(net[,1])))
{
Aindex=net[i,1]
Bindex=net[i,2]
counter=counter+1
		SIRENscore[counter]=0
		#SIRENscore[counter,1]=Aindex
		#SIRENscore[counter,2]=Bindex
		#SIRENscore[counter,4]=1
		prAB=PabCalculate(bMatrix[,,Aindex],bMatrix[,,Bindex])
		prA=prMatrix[Aindex,]

		prB=prMatrix[Bindex,]

		#Calculating Score for genes A and B
		#10: the number of bins for zero and non-zero values
		for (m in seq(10))
		{   
		    for (mm in seq(10))
	            {
			xN=log((prAB[m,mm])/(prA[m]*prB[mm]))
			tryCatch({
				if(xN>0)
		       		{
					 SIRENscore[counter]=SIRENscore[counter]+(prAB[m,mm])*weightMatrix[m,mm]*xN
				}
			 }
			 ,error=function(e) {}
			 ,finally={})
		    }
		}
}
}
return(SIRENscore)
}


SIREN=function(expMatrix,weightMatrix,net=NULL)
{
#this function in the main function of all the program
#exp=expression data

bMatrix=BMatrixFn(expMatrix)
prMatrix=PaCalculate(bMatrix)
SIRENmatrix=ScoreFn(bMatrix,weightMatrix,prMatrix,net)
#Pvalmatrix=SIRENpval(bMatrix,weightMatrix,SIREN)
cat("Preparing results...\n")
#result=as.data.frame(matrix(0,nrow=length(SIRENmatrix),ncol=4))
#names(result)=c('FirstNode','SecondNode','score','Pvalue')
result=as.data.frame(matrix(0,nrow=length(SIRENmatrix),ncol=3))
names(result)=c('FirstNode','SecondNode','score')

if(is.null(net))
{
counter=0
for(i in seq(length(bMatrix[1,1,])))
{
	for(j in seq(i+1,length(bMatrix[1,1,]),1))
	{
		counter=counter+1
		result[counter,1]=i
		result[counter,2]=j
		result[counter,3]=SIRENmatrix[counter]
		#result[counter,4]=Pvalmatrix[counter]
	}
}
}
else
{
for(i in seq(length(net[,1])))
{
	result[i,1]=net[i,1]
	result[i,2]=net[i,2]
	result[i,3]=SIRENmatrix[i]
	#result[counter,4]=Pvalmatrix[i]
	
}
}
return(result)
}


