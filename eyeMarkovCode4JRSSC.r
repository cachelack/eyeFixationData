
# Library for KNN classification
library(class)

# Library to display PNG images
library(png)

library(MASS);

BW_SCALE = 1.0;
GRID_SIZE= 100;
MCDRAWS  = 100;

# Change to wherever the data is kept on your system
#dataDir = "./";
dataDir = "/home/adam/Documents/RESEARCH/RSS2015/files-master/";

# Read in the entire data table
if(!exists("eyeData")){
  eyeData <- read.table(
    paste(dataDir,'data.csv',sep=""),h=TRUE,sep=','
  );
}

# some dimensions for ploting the images
eyeXlim <- c(-1280,1280)/2;
eyeYlim <- c(-1024,1024)/2;
eyeXlim.landscape <- c(-1024,1024)/2;
eyeYlim.landscape <- c(-768,768)/2;
eyeXlim.portrait <- c(-768,768)/2;
eyeYlim.portrait <- c(-1024,1024)/2;

# Returns subsets of data table

rss_getData <- function( imag=F, type=F, subj=F, fixt=F )
{
  if(type!=F)
  {
    chr = switch(
      type,
      normal   = "n",
      abnormal = "a",
      grayscale= "g"
    );
  }
  return( 
    subset( 
      eyeData, 
      tern(imag==F,T,image==imag) & 
      tern(type==F,T,condition==type) &
      tern(subj==F,T,subject==paste(chr,subj,sep="")) & 
      tern(fixt==F,T,fixation==fixt) 
    ) 
  );
}

# Ternary operator

tern <- function( a, b, c)
{
  if(a){ 
    return(b);
  } else {
    return(c);
  }
}

#
#  Add a cumulative duration column to the data table
#

if( ncol(eyeData)==8 ){
  dat = c();
  pb  = txtProgressBar(style=3,min=0,max=60);
  for( i in 1:60 ){
    for( j in c("normal","abnormal","grayscale") ){
      for( k in 1:10 ){
        tmp = rss_getData( i,j,k );
        if(nrow(tmp)==0) next;
        tmp = tmp[order(tmp[,4]),];
        tmp = cbind(tmp,cdur=c(cumsum(tmp[,3])));
        #tmp = cbind(tmp,cdur=c(0,cumsum(tmp[,3])[-nrow(tmp)]));
        dat = rbind(dat,tmp);
      }
    }
    setTxtProgressBar(pb,i);
  }
  close(pb);
  eyeData=dat; 
}

#########
#  Functions to fit the model from Section two
#  computing Bayes factors different cluster values
#########

#
#  Compute optimal clusters for all
#  img in images
#

rss_computeAllBayes <- function( 
  images, CLUSTTYPE=1, dst="eucl",meth="complete", kn=9
){
  x11(width=12,height=6);
  n  = length(images);
  # Initialize array of Bayes Factors
  # and array of clusters
  allBF = c();
  bf = matrix(0,n,3);
  cl = matrix(1,n,3);
  lk = matrix(0,n,30);
  offs  = min(images)-1;
  sTime = Sys.time();
  for( img in images )
  {
    print(paste("Testing Image ",img,"...",sep=""));
    nOut = rss_mmBayesFactor( img, 'normal',    PRINT=FALSE,
      CLUSTTYPE=CLUSTTYPE,dst=dst,meth=meth,kn=kn );
    aOut = rss_mmBayesFactor( img, 'abnormal',  PRINT=FALSE,
      CLUSTTYPE=CLUSTTYPE,dst=dst,meth=meth,kn=kn );
    gOut = rss_mmBayesFactor( img, 'grayscale', PRINT=FALSE,
      CLUSTTYPE=CLUSTTYPE,dst=dst,meth=meth,kn=kn );
    bf[img-offs,1] = min( nOut$bf );
    if(bf[img-offs,1]<.15) cl[img-offs,1] = which.min( nOut$bf );
    bf[img-offs,2] = min( aOut$bf );
    if(bf[img-offs,2]<.15) cl[img-offs,2] = which.min( aOut$bf );
    bf[img-offs,3] = min( gOut$bf );
    if(bf[img-offs,3]<.15) cl[img-offs,3] = which.min( gOut$bf );
    lk[img-offs,]  = c(
        nOut$like[nOut$kBest,],
        aOut$like[aOut$kBest,],
        gOut$like[gOut$kBest,]
    );
    print("Results: ");
    print(paste(" Bayes: ",bf[img-offs,]));
    print(paste(" Clust: ",cl[img-offs,]));
    allBF = rbind(allBF,nOut$bf,aOut$bf,gOut$bf);
  }
  print(Sys.time()-sTime);
  print( paste(
    "Average Time per pic:",(Sys.time()-sTime)/(3*n)
  ) );
  return( list(img=images,bf=bf,cl=cl,lk=lk,allBF=allBF) );
}

#
#  Compute best clusters for the Markov model 
#  by computing the Bayes factors for a 
#  specific image
#

rss_mmBayesFactor <- function( n,type, 
   DEBUG=FALSE, PRINT=TRUE, PLOT=TRUE, REMOVE=0,
   CLUSTTYPE=1,dst="euclidean",meth="complete",kn=3,
   cdr=c(0,5000), POSTHOC=FALSE
)
{
  # Get data
  dat   = rss_getData(n,type);
  if(ncol(dat)==9){
    dat = subset(dat,(cdr[1]<=dat[,9])&(dat[,9]<=cdr[2]));
  }
  print(paste("Data Size: ",nrow(dat)));
  for( s in REMOVE )
    dat = subset( dat, subject!=subj(type,s) );
  sCnt  = length(setdiff(1:10,REMOVE));
  # Max number of clusters to test
  maxK = 10;
  # matrix of likelihoods
  like1= matrix(-1000,maxK,10);
  # List of Markov posteriors
  mrkTran = list();
  mrkInit = list();
  for( k in 1:maxK ) # For each cluster count
  {
    mrkInit[[k]] = rep(0,k);
    mrkTran[[k]] = matrix(0,k,k);
    if(PRINT)
      print(paste("Testing k-means ",k,"...",sep=''));
    for( s in setdiff(1:10,REMOVE) ) # For each subject (cross validate)
    {
      # Get Training and Test Data 
      train = subset( dat, subject!=subj(type,s) );
      test  = subset( dat, subject==subj(type,s) );
      if( nrow(test)<2 ) next;
      test  = rss_sortFix(test);
      # If no clustering (i.e. clusters==1 )
      if( k==1 ){
        dTotal = kde2d(
          train$fx,train$fy,lims = c(eyeXlim, eyeYlim),n=GRID_SIZE
        );
        like1[k,s] = calcLike0( dTotal, test );
        next;
      }
      # Run k means
      kdat  = cbind( train$fx, train$fy );
      colnames(kdat) <- c("x","y");
      if(CLUSTTYPE==1){
        km    = kmeans(kdat,k,iter.max=20,nstart=10);
      } else if( CLUSTTYPE==2 ){
        hc = hclust( dist(kdat,method=dst),method=meth );
        ct = as.double(cutree( hc, k=k ));
        km = list(cluster=ct,centers=c(),size=c());
        for( kk in 1:k){
          km$size = c(km$size,length(which(ct==kk)));
          km$centers = rbind(km$centers,
            apply( subset(kdat,ct==kk),2,mean )
          );
        }
      }
      # Return if bad clusters
      #print(paste(km$size,sum(km$size)));
      if( sum( km$size<3 )>0 ) next;
      # Assign cluster to Train and Test sets
      train = cbind( train, cluster=km$cluster );
      test  = cbind( test,  cluster=0 );
      # Count number of transitions between clusters
      # Use to update Dirichlet Prior
      markv = rss_constrMarkov( n,k,type,train,km );
      mrkInit[[k]] = mrkInit[[k]] + markv$iProb;
      mrkTran[[k]] = mrkTran[[k]] + markv$tProb;
      # Relabel Clusters
      km$cluster = markv$km$cluster;
      km$centers = markv$km$centers;
      train$cluster = km$cluster;
      # Pick best cluster for each element of test set
      #for( i in 1:nrow(test) )
      #  test$cluster[i] = chooseClust( test$fx[i],test$fy[i],train );
      test$cluster = chooseClust2( test, train, kn );
      # Construct empirical density estimate for each cluster
      dClust = constrClustDen( train, k );
      # Compute likelihood given transition probabilities and
      # density estimates
      like1[k,s] = calcLike1( dClust, markv, km, test, MCDRAWS );
      if(DEBUG)
      {
        print(paste("Like: ",like1[k,s]));
        rss_plotTrainTest(n,k,type,train,test,km);
        readline();
      }
    }
    #mrkInit[[k]] = mrkInit[[k]]/10;
    #mrkTran[[k]] = mrkTran[[k]]/10;
  }
  
  # Choose best number of clusters excluding k=1
  kBest = which.max( rowSums((like1[-1,])) )+1;
  # Average over all subjects
  kLike = (rowSums((like1))/sCnt);
  # Compare k>1 clusters to k=1 (i.e. Bayes Factor)
  #bf    = pmin(exp(kLike[1]-kLike),rep(1));
  bf    = exp(kLike[1]-kLike);
  if( min(bf>0.15) ) kBest=1;
  if(PLOT)
  {
    if(kBest==1){
      rss_plotImageData(n,type);
    } else {
      clDat = rss_clusterDensity(n,kBest,type,min(bf[kBest]),PRINT);
    }
  }
  if(PRINT){
    print(signif(kLike,4));
    print(signif(pmin(exp(kLike[1]-kLike),rep(1)),4));
    print(signif(like1[kBest,],4));
  }

  if(POSTHOC)
    return( posthocMarkov(n,type,kBest,2000,JUSTCNT=TRUE) );

  return(list(
    like=like1,bf=bf,kBest=kBest,
    mInit=mrkInit[[kBest]],#/sum(mrkInit[[kBest]]),
    mTran=mrkTran[[kBest]]#/rowSums(mrkTran[[kBest]])
  ));
}

#
# Compute Likelihood in null case
# that is, the random draw model without
# any Markov transitions.
#

calcLike0 <- function( den, dat )
{
  n = nrow(dat);
  like = 0
  for( i in 1:n )
  {
    # Find matrix coord cooresponding to actual coord
    x = max(1,which(den$x < dat$fx[i]));
    y = max(1,which(den$y < dat$fy[i]));
    like = like + log(den$z[x,y]);
  } 
  return(like);
}

#
# Compute Likelihood in given kmean clusters km and
# posterior hyperparameters for the Markov chain.
#

calcLike1 <- function( dens, markv, km, dat, mcDraws )
{
  n  = nrow(dat);
  likeTot = 0;
  # MC integration of Bayes factor
  # over the Dirichlet distribution
  for( m in 1:mcDraws )
  {
    # Draw initial and transition probabilities
    # from the Dirichlet dist
    iProb = rdirich( markv$iProb );
    tProb = rdirich( markv$tProb );
    # Calc initial fixation probability
    i  = 1;
    cl = dat$cluster[i];
    x  = max(1,which(dens[[cl]]$x < dat$fx[i]));
    y  = max(1,which(dens[[cl]]$y < dat$fy[i]));
    like = log(iProb[cl]) + log(dens[[cl]]$z[x,y]);
    for( i in 2:n ) # for each subsequent fixation
    {
      # Calc fixation likelihood and transition probability
      old= dat$cluster[i-1];
      cl = dat$cluster[i];
      x  = max(1,which(dens[[cl]]$x < dat$fx[i]));
      y  = max(1,which(dens[[cl]]$y < dat$fy[i]));
      like = like + log(tProb[old,cl]) + log(dens[[cl]]$z[x,y]);
    }
    likeTot = likeTot + exp(like);
  } 
  return(log(likeTot)-log(mcDraws));
}

#
#  Construct a list of density estimates,
#  one for each cluster
#

constrClustDen <- function(train,k)
{
  dens  = list();
  for( i in 1:k )
  {
    clust     = subset(train,cluster==i);
    #print(paste(bandwidth.nrd(clust$fx),bandwidth.nrd(clust$fy)));
    #print(rbind(clust$fx,clust$fy));
    dens[[i]] = kde2d(
      clust$fx,clust$fy,lims = c(eyeXlim, eyeYlim),n=GRID_SIZE
    );
  }
  return(dens);
}

#
#  Use observed data to update Dirichlet priors
#

rss_constrMarkov <- function( n,k,type,dat,km )
{
  chr = switch(
    type,
    normal   = "n",
    abnormal = "a",
    grayscale= "g"
  );

  trial = list();
  # Jeffreys Priors
  initP = rep(1/2,k);
  trans = matrix(1/2,k,k);
  for( i in 1:10 )
  {
    # Look at ith subject
    trial[[i]] = subset(dat,dat$subject==paste(chr,i,sep=""));
    trial[[i]] = rss_sortFix(trial[[i]]);
    # Skip if no data
    if(nrow(trial[[i]])==0) next;
    # Update Priors
    initP[trial[[i]]$cluster[1]] = initP[trial[[i]]$cluster[1]]+1;
    for( j in 2:nrow(trial[[i]]) )
      trans[trial[[i]]$cluster[j-1],trial[[i]]$cluster[j]] =
        trans[trial[[i]]$cluster[j-1],trial[[i]]$cluster[j]]+1;
  }
  
  # Order probabilities and clusters based
  # on the most populous cluster
  #ord   = order(km$size,decreasing=TRUE);
  ord   = order(km$centers[,1],decreasing=FALSE);
  ord2  = order(ord);
  initP = initP[ord];
  trans = trans[ord,];
  trans = trans[,ord];
  km$cluster = ord2[km$cluster];  
  km$centers = km$centers[ord,];

  #initP = initP/sum(initP);
  #trans = trans/rowSums(trans);

  return(list(iProb=initP,tProb=trans,km=km)); 
}

#
#  Plotting Functions
#

rss_plotCluster <- function(n,k,type,kdat,km,NOTITLE=FALSE)
{
  img = rss_readImage(n,type);
  mycol = rainbow(k);
  if(NOTITLE){
    plot( 
      0, 0, xlim=eyeXlim*5/6, ylim=eyeYlim*4/5, xlab="",ylab="",
      type="n",xaxt='n',yaxt='n'
    );
  } else {
    plot( 
      0, 0, xlim=eyeXlim*5/6, ylim=eyeYlim*4/5, xlab="",ylab="",
      main=paste("Fixations on image",n),type="n",
      xaxt='n',yaxt='n'
    );
  }
  rasterImage( img , 
    eyeXlim.landscape[1], eyeYlim.landscape[1], 
    eyeXlim.landscape[2], eyeYlim.landscape[2]
  );
  points( kdat, pch=19, col=mycol[km$cluster], cex=1.1 );
  points( kdat, pch=1,  col="black", cex=1.1 );
  points( km$centers,   col=mycol, pch = 13, cex = 2)
}

rss_plotTrainTest <- function(n,k,type,train,test,km)
{
  img = rss_readImage(n,type);
  mycol = rainbow(k);
  plot( 
    0, 0, xlim=eyeXlim, ylim=eyeYlim, xlab="x",ylab="y",
    main=paste("Fixations on image",n),type="n"
  );
  rasterImage( img , 
    eyeXlim.landscape[1], eyeYlim.landscape[1], 
    eyeXlim.landscape[2], eyeYlim.landscape[2]
  );
  points( train$fx,train$fy, col=mycol[train$cluster], cex=1.1 );
  points( test$fx,test$fy, pch=19, col=mycol[test$cluster], cex=1.1 );
  points( test$fx,test$fy, pch=1,  col="black", cex=1.1 );
  points( km$centers,   col=mycol, pch = 13, cex = 2)
}


rss_clusterDensity <- function(n,k,type,bf,PRINT=TRUE)
{
  dat = rss_getData(n,type);
  kdat= cbind(dat$fx,dat$fy);
  colnames(kdat) <- c("x","y");
  km = kmeans(kdat,k,iter.max=20,nstart=10);
  ## reorder
  ord   = order(km$centers[,1],decreasing=FALSE);
  ord2  = order(ord);
  km$cluster = ord2[km$cluster];  
  km$centers = km$centers[ord,];
  ##
  print(paste("Type:",type));
  if(PRINT){
    for( i in 1:k )
      print(paste("  Cluster ",i,": ",sum(km$cluster==i),sep=""));
  }

  dat  = cbind(dat,cluster=km$cluster);
  dens = constrClustDen(dat,k);
  rss_plotClusterDensity(n,k,type,dat,km,dens);
  title(
    paste("Fixations on Image ",n,", Bayes Factor = ",signif(bf,3),sep=""),
    line=-2.5,outer=TRUE,cex.main=3
  );
  return(dat);
}

rss_plotClusterDensity <- function(n,k,type,dat,km,dens)
{
  mm  = length(dens);
  if(mm==1){
    par(mfrow=c(1,2));
  } else if(mm==2) {
    layout(matrix(c(1,1,1,1,2,3), 2, 3, byrow = FALSE));
  } else if(mm<5) {
    layout(matrix(c(1,1,1,1,2,4,3,5), 2, 4, byrow = FALSE));
  } else if(mm<7) {
    layout(matrix(c(1,1,1,1,2,5,3,6,4,7), 2, 5, byrow = FALSE));
  } else {
    layout(matrix(c(rep(1,9),2,5,8,3,6,9,4,7,10), 3, 6, byrow = FALSE));
  }
  rss_plotCluster(n,k,type,dat,km,NOTITLE=TRUE);
  for( i in 1:k )
  {
     clust     = subset(dat,cluster==i);
     dens[[i]] = kde2d(
       clust$fx,clust$fy,lims = c(eyeXlim, eyeYlim),n=100
     );
  }
  for( i in 1:k )
  {
     image(
       dens[[i]]$x,dens[[i]]$y,log(dens[[i]]$z),
       zlim=max(log(dens[[i]]$z))+c(-10,0),
       xlab='',ylab='',xaxt='n',yaxt='n'
     );
     for( j in 1:k )
       image(dens[[j]]$x,dens[[j]]$y,log(dens[[j]]$z),
         zlim=max(log(dens[[j]]$z))+c(-10,0),
         add=TRUE,col=cm.colors(12,alpha=0.5)
       );
     image(
       dens[[i]]$x,dens[[i]]$y,log(dens[[i]]$z),
       zlim=max(log(dens[[i]]$z))+c(-10,0),add=TRUE
     );
  }
}



##########
#  Functions related to reading in the data.
##########

#
# Read in Image of #n of type 'type'
# where 'type' is one of "normal", "grayscale", "abnormal" 
#

rss_readImage <- function( n, type )
{
  typeCode = switch(
    type,
    normal   = "c",
    abnormal = "ab",
    grayscale= "bw"
  );
  return( 
    readPNG(
      paste(
        dataDir,'images/',type,'/',n,'_',typeCode,".png",sep=""
      ) 
    )
  );
}

#
#  Extract data subset of 
#  image     imag in {1,...,60} 
#  condition type in {"normal","grayscale","abnormal"}
#  subject   subj in {1,...,10}
#  fixation  fixt in {1,...,26}
#

rss_getData <- function( imag=F, type=F, subj=F, fixt=F )
{
  if(type!=F)
  {
    chr = switch(
      type,
      normal   = "n",
      abnormal = "a",
      grayscale= "g"
    );
  }
  return( 
    subset( 
      eyeData, 
      tern(imag==F,T,image==imag) & 
      tern(type==F,T,condition==type) &
      tern(subj==F,T,subject==paste(chr,subj,sep="")) & 
      tern(fixt==F,T,fixation==fixt) 
    ) 
  );
}

#
#  Computes Euclidean distance between adjacent fixations
#

rss_computeSacs <- function()
{
  dat = rss_getData();
  dat = dat[ order(dat$subject,dat$image,dat$fixation), ];
  n   = nrow(dat);
  sac = sqrt( 
    (dat$fx[2:n]-dat$fx[1:(n-1)])^2 + 
    (dat$fy[2:n]-dat$fy[1:(n-1)])^2 
  );
  indx= dat$image[1:(n-1)]==dat$image[2:n];
  sac = c(0,sac);
  indx= c(0,indx);
  dat$sac = sac*indx;
 
  return(dat);
}

#
#  Similar to rss_getData, but returns Euclidean 
#  distance between successive fixations
#

rss_getSacData <- function( imag=F, type=F, subj=F, fixt=F )
{
  sacData = rss_computeSacs();
  if(type!=F)
  {
    chr = switch(
      type,
      normal   = "n",
      abnormal = "a",
      grayscale= "g"
    );
  }
  return( 
    subset( 
      sacData, 
      tern(imag==F,T,image==imag) & 
      tern(type==F,T,condition==type) &
      tern(subj==F,T,subject==paste(chr,subj,sep="")) & 
      tern(fixt==F,T,fixation==fixt) 
    ) 
  );
}
#
# Sort by fixation
#

rss_sortFix <- function( dat )
{
  return( dat[order(dat$fixation),] );
}

#
# Read in image and plot data
#

rss_plotImageData <- function( n, type  )
{
  img = rss_readImage(n,type);
  dat = rss_getData(imag=n,type=type);
  plot( 
    0, 0, xlim=eyeXlim, ylim=eyeYlim, xlab="x",ylab="y",
    main=paste("Fixations on image",n),type="n"
  );
  rasterImage( img , 
    eyeXlim.landscape[1], eyeYlim.landscape[1], 
    eyeXlim.landscape[2], eyeYlim.landscape[2]
  );
  points( dat$fx, dat$fy, pch=19, col="red", cex=1.1 );
}

rss_plotMultiPics <- function( images, type )
{
  n = length(images);
  d = ceiling(sqrt(n));
  par(mfrow=c(d,d));
  for( i in images )
    rss_plotImageData( i, type );
}

#
# Misc
#

# Subject String

subj <- function(type,s)
{
  chr=switch(
    type,
    normal   = "n",
    abnormal = "a",
    grayscale= "g"
  );
  return(paste(chr,s,sep=""));
}

# Ternary operator

tern <- function( a, b, c)
{
  if(a){ 
    return(b);
  } else {
    return(c);
  }
}

#
# Choose best Cluster
#

# use k-nearest-neighbours
chooseClust2 <- function( test, train, k=3 )
{
  res = knn( 
    cbind(train$fx,train$fy),cbind(test$fx,test$fy), 
    train$cluster,k=k
  );
}

# use distance to center
chooseClust <- function( x, y, train )
{
  dst = rep(0,nrow(km$centers));
  for( i in 1:nrow(km$centers) )
  {
    dst[i] = 
      (x-km$centers[i,1])^2 + (y-km$centers[i,2])^2;
  }
  return(which.min(dst));
}

#
#  Random Dirichlet draw
#

rdirich <- function( alf )
{
  n   = length(alf)
  out = rgamma(n,alf,1);
  if(is.matrix(alf)){
    out= matrix(out,nrow(alf),ncol(alf));
    return( out/rowSums(out) ); 
  } else {
    return( out/sum(out) );
  }
}

#
#  Row/Column min/max
#

rMin <- function(mat){
  return(do.call(pmin, as.data.frame(mat)));
}

rMax <- function(mat){
  return(do.call(pmax, as.data.frame(mat)));
}

cMax <- function(mat){
  return(do.call(pmax, as.data.frame(t(mat))));
}

cMin <- function(mat){
  return(do.call(pmin, as.data.frame(t(mat))));
}

######################################
#### Extra
######################################

rss_plotImageData2 <- function( n, type, NOAXIS=TRUE  )
{
  img = rss_readImage(n,type);
  dat = rss_getData(imag=n,type=type);
  if(NOAXIS){
    plot( 
      0, 0, xlim=eyeXlim*5/6, ylim=eyeYlim*4/5, xaxt='n',yaxt='n',
      xlab="",ylab="",type="n"
    );
  } else {
    plot( 
      0, 0, xlim=eyeXlim*5/6, ylim=eyeYlim*4/5, xlab="x",ylab="y",
      type="n"
    );
  }
  rasterImage( img , 
    eyeXlim.landscape[1], eyeYlim.landscape[1], 
    eyeXlim.landscape[2], eyeYlim.landscape[2]
  );
  points( dat$fx, dat$fy, pch=19, col="red", cex=1 );
  points( dat$fx, dat$fy, pch=1, col="black", cex=1 );
}

rss_plotImageClusters2 <- function( 
  n, type, k, meth="complete", dst="euclidean" 
){
  img = rss_readImage(n,type);
  dat = rss_getData(imag=n,type=type);
  
  # k-means
  kdat = cbind(dat$fx,dat$fy);
  colnames(kdat) = c("x","y");
  km = kmeans( kdat,k,iter.max=20,nstart=10 );
  dat$kClust = km$cluster;
  
  # hierachical
  hc = hclust( dist(cbind(dat$fx,dat$fy),method=dst),method=meth );
  dat$hClust = cutree(hc, k = k);

  par(mfcol=c(1,2));
  plot( 
    0, 0, xlim=eyeXlim*5/6, ylim=eyeYlim*4/5, xlab="x",ylab="y",
    type="n"
  );
  rasterImage( img , 
    eyeXlim.landscape[1], eyeYlim.landscape[1], 
    eyeXlim.landscape[2], eyeYlim.landscape[2]
  );
  points( dat$fx, dat$fy, pch=19, col=dat$kClust+1, cex=1 );
  points( dat$fx, dat$fy, pch=1, col="black", cex=1 );
  title("k-means Clustering",cex.main=2);

  plot( 
    0, 0, xlim=eyeXlim*5/6, ylim=eyeYlim*4/5, xlab="x",ylab="y",
    type="n"
  );
  rasterImage( img , 
    eyeXlim.landscape[1], eyeYlim.landscape[1], 
    eyeXlim.landscape[2], eyeYlim.landscape[2]
  );
  points( dat$fx, dat$fy, pch=19, col=dat$hClust+1, cex=1 );
  points( dat$fx, dat$fy, pch=1, col="black", cex=1 );
  title("Hierarchical Clustering",cex.main=2);
}

