#This version May 23, 2023

#### arguments from bash file
args=commandArgs(trailingOnly = TRUE)
#1: model
#2: pedigree filename (optional for GBLUP)
#3: vcf filename
#4: phenotype filename
#5: allele frequencies (optional)

##### Input files and model specification

#model
model=args[1]#accepts GBLUP, ssGBLUP
print(paste("Fitting", model))

#file name of pedigree
pedigreename=args[2]
#contains an object with 3 columns: id, father id, mother id
if(exists("pedigreename",envir=.GlobalEnv)){
  ped0<-read.table(pedigreename,header=TRUE,sep=",")
  print(paste("Pedigree file:", pedigreename))
}

#file name of vcf
vcfname=args[3]
print(paste("vcf file:", vcfname))

#file name of environmental effects, if any
environame="simenvirofull"

allelefreq=args[5]
if(allelefreq!='NA'){ #is there a given name of existing allele frequencies?
  load(allelefreq)
  #should contain a vector called 'allelefreqknown'
}


#filename of phenotypes
phenoname=args[4]
print(paste("Phenotype file:", phenoname))

dataobs<-read.table(phenoname,header=TRUE,sep=",")[,-1]
Y<-dataobs[,2]
if(dim(dataobs)[2]>2){ #extra columns past tags, pheno, are enviro covariates
  X=dataobs[,-(1:2)]
  if(any(X[,1]!=1)){
   X<-cbind(rep(1,length(Y)),X)
  }
  fe=dim(X)[2]
} else{
  X=matrix(data=1,nrow=length(Y))
  fe=1
}

#option to provide heritability if known. Otherwise estimated as part of model
h2known=FALSE
if(h2known==TRUE){
  h2=0.5
}


#####Converting vcf file

library(vcfR)

print("Reading vcf file...")

vcf <- read.vcfR(file=vcfname, verbose = FALSE)
variants<-dim(vcf)[1]

print("DONE")

idtags<-dimnames(vcf@gt)[[2]][-1]
genotags<-idtags
if(any(duplicated(idtags))){
  stop("Multiple records given for same individual")
}

#vcf<-as.data.frame(vcf@gt)
vcf<-as.matrix(vcf@gt,nrow(variants))
vcf<-vcf[,-1]

#####when testing partial geno:

print("Converting vcf file...")
convfunc <- function(x){ #x is a vector
  temp=NA
  for(i in 1:length(x)){
    if(x[i]=="0|0"){temp[i]=0}
    else if(x[i]=="0|1"){temp[i]=1}
    else if(x[i]=="1|0"){temp[i]=1}
    else if(x[i]=="1|1"){temp[i]=2}
    else{temp[i]=-1} #if have misspecified or multi-allelic
  }
  return(temp)
}

genomat<-matrix(unlist(lapply(vcf,convfunc)),nrow=variants)
colnames(genomat)<-genotags
print("DONE")

print("Calculating relationship matrix...")
#function to calculate E[M] for each locus, ignoring missing values
Mexpfunc<- function(x){ #x is a vector
  return(mean(x[x>-1]))
}

Mvarfunc<- function(x){ #x is a vector
  return(var(x[x>-1]))
}


EM<-apply(genomat,1,Mexpfunc)
VM<-apply(genomat,1,Mvarfunc)
#only need loci where there is actually some variation, i.e. VM>0
stor<-which(EM<0.05|EM>1.95|is.na(EM)|VM==0) #loci to be removed
genomat<-genomat[-stor,]
EM<-EM[-stor]
VM<-VM[-stor]

#if allele frequencies are provided use those:
if(exists("allelefreqknown")){
  stop("to do")
  EM=allelefreqknown[-stor]
  #still have same locus filtering as above (based on observed alleles)
}


#calculate 2pq for each locus
pq2=EM*(1-EM/2)

###documentation: only keep loci where minor allele is at least 5%


#set missing values to 0 for that locus and adjust others by EM
for(i in 1:dim(genomat)[1]){
  temp=which(genomat[i,]==-1)
  if(length(temp)>0){
    genomat[i,temp]<-0
    
    genomat[i,-temp]<-genomat[i,-temp]-EM[i]
  } else {
    genomat[i,]<-genomat[i,]-EM[i]
  }
  
  
}
#change this so not a for loop? 'apply' with some function


#creation of incidence matrix Z

#presently, assume that all phenotyped individuals have known environmental covariates

ngen=dim(genomat)[2]
nind=length(unique(c(dataobs[,1],names(vcf[1,]))))

Z=matrix(data=0,nrow=length(Y),ncol=nind)


nobs=length(unique(dataobs[,1]))

if(length(Y)==dim(genomat)[2]){
  if(model=='ssGBLUP'){ #warning that pedigree info unused when all individuals are genotyped, use GBLUP
    warning("All individuals are genotyped and pedigree information is not required. Fitting GBLUP instead of ssGBLUP.")
    model='GBLUP'
    
    
  }

  
  
} else if(length(Y)<dim(genomat)[2]){
  
} else if (model=='GBLUP'){ #throw warning if ungenotyped individuals in GBLUP: only predict genotyped ind
  stop("GBLUP requires that all individuals are genotyped.")
  #warnings("GBLUP requires all individuals are genotyped. Prediction will only use genotyped individuals.")
}

for(ind in 1:length(Y)){
  name=dataobs[,1][ind]
  col=which(idtags==name)
  
  #only execute following line if dataobs[,1][ind] is actually in vcf names (as idtag may have names from ped, not geno)
  Z[ind,col]=1
}

##### check if first column of X is intercept, add intercept if not
if(!all(X[,1]==1)){
  X<-cbind(rep(1,length(Y)),X)
}

    

###create A matrix if ssGBLUP or if GBLUP with pedigree but no allelefrequencies

    
    
if(model=='ssGBLUP'|all(model=='GBLUP'|!exists("allelefreqknown")|exists("ped"))){
  library("pedigreemm")
  
  #ped0 is expected to be idtags, but 'pedigree' command needs numbers
  #convert from idtag to position in unique(idtag)
  idtags<-unique(c(idtags,ped0[,1]))
  pedtags<-ped0[,1]
  
  
  
  ped<-pedigree(sire=ped0[,2],dam=ped0[,3],label=ped0[,1])
  A<-getA(ped)
  Ainv<-getAInv(ped)
  
  
  #use a vector for this? A[idtag,idtag] or similar
  #assume dim(A) >= dim(G) ? i.e. no geno ind missing in ped
  #add them to ped if so
  
  
}





##### GBLUP and ssGBLUP code
library(rrBLUP)

#create G matrix
G=t(genomat)%*%(genomat)/sum(pq2) #will be slow doing it this way with many ind.. is there a faster option?



##### GBLUP code
if(model=='GBLUP'){
#use 'genotags' sort order as all ind are genotyped

if(!exists("allelefreqknown")){
  if(exists("A")){ #if have known A matrix, use 0.99G+0.01A
    G=0.99*G+0.01*A[genotags,genotags]
  } else if(min(eigen(G)$values)<=0){
    stop("Matrix inversion issues as allele frequencies are estimated from observed genotypes. To fix this, provide either existing allele frequencies (preferred) or pedigree information.")
  } else{
    warning("Matrix inversion issues may have affected prediction as allele frequencies are estimated from observed genotypes. To fix this, provide either existing allele frequencies (preferred) or pedigree information.")
  }
  
}

print("DONE")
print("Fitting GBLUP")

GBLUP <- mixed.solve(y=Y, Z=Z, K=as.matrix(G), X=X, method="REML",
                  bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
print("DONE")
predictions<-GBLUP$u
rownames(predictions)<-genotags

}


##### ssGBLUP code
if(model=='ssGBLUP'){
#use 'pedtags' sort order as all ind are in pedigree
  
#determine vector of genotyped individuals: which ind in pedtags' are in 'genotags' ?
geno=rep(0,nind)
for(i in 1:nind){
  if(any(genotags==pedtags[i])){
    geno[i]=1
  }
}
  
  
G=0.99*G+0.01*A[genotags,genotags]

Ginv=solve(G)

Hinv=Ainv
A22inv=solve(A[genotags,genotags])

Hinv[genotags,genotags]<-Hinv[genotags,genotags]+Ginv-A22inv
H=solve(Hinv)

print("DONE")
print("Fitting ssGBLUP")
ssGBLUP <- mixed.solve(y=Y, Z=Z, K=as.matrix(H), X=X, method="REML",
                     bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
print("DONE")
    
predictions<-ssGBLUP$u
rownames(predictions)<-pedtags

}

##### output
print("Predicted values:")
predictions