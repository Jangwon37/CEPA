library(survival);
library(MLEcens);
source('../cOPT/R/opt2d.R');
source('../cOPT/R/dabrowska.R');
source('../cOPT/R/linying.R');
source('./Script.JW.01.R');

# functions
rowVars = function(E,na.rm=FALSE) {
  n = ncol(E);
  m = rowMeans(E,na.rm=na.rm);
  mm = rowMeans(E^2,na.rm=na.rm);
  v = (n*mm - n*m^2)/(n-1);
  return( v );
}

random_data_surv = function(N,lt=1,lc=1/2,of=0) {
  # T1 = rexp(N,lt); T2 = T1 + rexp(N,lt) + rnorm(N,0,1); T3 = T2 + rexp(N,lt) + rnorm(N,0,1);
  T1 = rexp(N,lt); T2 = T1 + rexp(N,lt); T3 = T2 + rexp(N,lt);
  C1 = rexp(N,lc)+of; C2 = rexp(N,lc)+of; C3 = rexp(N,lc)+of
  X1 = Surv(apply(cbind(T1,C1),1,min),(T1<=C1) );
  X2 = Surv(apply(cbind(T2,C2),1,min),(T2<=C2) );
  X3 = Surv(apply(cbind(T3,C3),1,min),(T3<=C3) );
  X = cbind(as.data.frame(X1),X2, X3); colnames(X) = c('X1','X2','X3');
  A = c(0,max(X1[,1])*1.1,max(X2[,1])*1.1,max(X3[,1])*1.1);
  # X: bivariate survival times
  # A: region of interest
  # T: true event times
  return( list(X=X, A=A, T=cbind(T1,T2,T3)) );
}

folder_= '../../result/simulation_data/exp'
ifelse(dir.exists(folder_), 'Folder exists already', dir.create(folder_))

if (!file.exists(sprintf("%s/prAB_result.csv", folder_))) {
  prAB_result = data.frame()
  prAB_result_db = data.frame()
  prAB_result_ly = data.frame()
  prAB_result_mle = data.frame()
}

file_names <- list.files(path = folder_)
numeric_files <- as.numeric(gsub(pattern = "\\..*$", replacement = "", file_names))
numeric_files <- numeric_files[!is.na(numeric_files)]

if(length(numeric_files) == 0) {
  loop_start <- 1
} else if(length(numeric_files) == 100) {
  stop("finish")
} else {
  loop_start <- max(numeric_files)
}

for(loop in loop_start:100){
  
  ifelse(dir.exists(sprintf("%s/%s/",folder_,loop)), 'Folder exists already', dir.create(sprintf("%s/%s/",folder_,loop)))
  
  print(loop)
  
  simul_data_X = data.frame(random_data_surv(500)[1])
  simul_data_A = data.frame(random_data_surv(500)[2])
  simul_data_T = data.frame(random_data_surv(500)[3])
  
  outputfile_X = sprintf("%s/%s/output.random_data_sim_X.csv", folder_, loop);
  write.csv(simul_data_X, file=outputfile_X)
  
  outputfile_A = sprintf("%s/%s/output.random_data_sim_A.csv", folder_, loop);
  write.csv(simul_data_A, file=outputfile_A)
  
  outputfile_T = sprintf("%s/%s/output.random_data_sim_T.csv", folder_, loop);
  write.csv(simul_data_T, file=outputfile_T)
  
  RECOV.EVENT <- read.csv(sprintf('%s/%s/output.random_data_sim_X.csv', folder_, loop), row.names = 1)
  
  prAB_values = data.frame()
  prAB_values_db = data.frame()
  prAB_values_ly = data.frame()
  prAB_values_mle = data.frame()
  
  for(i in 1:2){
    for(j in (i+1):3){
      
      inputfile = sprintf("%s/%s/input.surv%s%s.txt",folder_,loop,i,j);
      
      # copt
      outputfile = sprintf("%s/%s/output.surv%s%s.txt", folder_,loop,i,j);
      
      DATA.SURV = cbind(RECOV.EVENT[,2*i-1], RECOV.EVENT[,2*i],
                        RECOV.EVENT[,2*j-1], RECOV.EVENT[,2*j])
      
      write.table(DATA.SURV, file = inputfile, sep='\t',
                  col.names = FALSE, row.names = FALSE)
      
      system(sprintf("../cOPT/copt ./%s -o ./%s",inputfile,outputfile));
      result.opt = read.table(outputfile, sep = '\t');
      
      prAB_value = PrAB(result.opt)[1]
      
      MIN = floor(min(simul_data_A));
      MAX = ceiling(max(simul_data_A));
      A = seq(from = MIN, to = MAX, by = (MAX-MIN)/5);
      X_surv = cbind(as.data.frame(Surv(RECOV.EVENT[,2*i-1], RECOV.EVENT[,2*i])), Surv(RECOV.EVENT[,2*j-1], RECOV.EVENT[,2*j]))
      
      # dabrowska
      
      outputfile_db = sprintf("%s/%s/dabrowska_output.surv%s%s.txt", folder_,loop,i,j);
      
      TM = NULL;
      for( aaa in A ) { for( bbb in A )  TM = rbind(TM,c(aaa,bbb)); }
      S = dabrowska(RECOV.EVENT[,2*i-1], RECOV.EVENT[,2*j-1], RECOV.EVENT[,2*i], RECOV.EVENT[,2*j], TM)
      ss = t(rbind(cbind(matrix(S[,3],nrow=length(A)),0),0));
      tlist1 = c(A,Inf);
      tlist2 = c(A,Inf);
      SS = NULL;
      for( ii in 1:(length(tlist1)-1) ) {
        for( jj in 1:(length(tlist2)-1) ) {
          p = ss[ii,jj] - ss[ii+1,jj] - ss[ii,jj+1] + ss[ii+1,jj+1];
          SS = rbind(SS,c(tlist1[ii],tlist1[ii+1],tlist2[jj],tlist2[jj+1],0,p));
        }
      }
      SS[,5] = SS[,6]/(SS[,2]-SS[,1])/(SS[,4]-SS[,3]);
      Sdb = SS;
      prAB_value_db = PrAB(Sdb,MAX,MAX)[1]
      
      if(prAB_value_db>1){prAB_value_db = 0.99}
      
      if(prAB_value_db<0){prAB_value_db = 0.01}
      
      # linying estimates
      
      outputfile_ly = sprintf("%s/%s/linying_output.surv%s%s.txt", folder_,loop,i,j);
      tlist1 = c(A,Inf);
      tlist2 = c(A,Inf);
      SS = NULL;
      for( ii in 2:(length(tlist1)) ) {
        for( jj in 2:(length(tlist2)) ) {
          p = linying(X_surv, tlist1[ii-1],tlist2[jj-1]) - linying(X_surv,tlist1[ii-1],tlist2[jj]) - linying(X_surv,tlist1[ii],tlist2[jj-1]) + linying(X_surv,tlist1[ii],tlist2[jj]);
          SS = rbind(SS,c(tlist1[ii-1],tlist1[ii],tlist2[jj-1],tlist2[jj],0,p));
        }
      }
      SS[,5] = SS[,6]/(SS[,2]-SS[,1])/(SS[,4]-SS[,3]);
      Sly = SS;
      prAB_value_ly = PrAB(Sly,MAX,MAX)[1]
      
      if(prAB_value_ly>1){prAB_value_ly = 0.99}
      
      if(prAB_value_ly<0){prAB_value_ly = 0.01}
      
      
      # MLEcens esimates
      epsilon = 1e-6
      outputfile_mle = sprintf("%s/%s/mlecens_output.surv%s%s.txt", folder_,loop,i,j);
      TTE = cbind(round(RECOV.EVENT[,2*i-1],digits = 2),round(RECOV.EVENT[,2*i], digits = 2),round(RECOV.EVENT[,2*j-1], digits = 2),round(RECOV.EVENT[,2*j], digits = 2));
      
      RR <- matrix(ncol = ncol(TTE), nrow = nrow(TTE))
      RR[,1] <- TTE[,1] * (1 - TTE[,2])
      RR[,2] <- TTE[,1] + TTE[,2]+ epsilon
      RR[,3] <- TTE[,3] * (1 - TTE[,4])
      RR[,4] <- TTE[,3] + TTE[,4]+ epsilon
      
      mle = computeMLE(R=RR, B=c(1,1,1,1));
      SS = NULL;
      for( ii in 1:nrow(mle$rects)){
        SS = rbind(SS,c(mle$rects[ii,1],mle$rects[ii,2],mle$rects[ii,3],mle$rects[ii,4],0,mle$p[ii]));
      }
      SS[,5] = SS[,6]/(SS[,2]-SS[,1])/(SS[,4]-SS[,3]);
      Smle = SS;
      
      prAB_value_mle = PrAB(Smle,MAX,MAX)[1]
      
      if (is.na(prAB_value_mle)) {
        prAB_value_mle <- 0.01
      } else {
        
        if(prAB_value_mle>1){prAB_value_mle = 0.99}
        if(prAB_value_mle<0){prAB_value_mle = 0.01}
      }
      
      
      if (i == 1 & j == 2) {
        prAB_values <- prAB_value
        prAB_values_db <- prAB_value_db
        prAB_values_ly <- prAB_value_ly
        prAB_values_mle <- prAB_value_mle
      } else {
        prAB_values <- cbind(prAB_values, prAB_value)
        prAB_values_db <- cbind(prAB_values_db, prAB_value_db)
        prAB_values_ly <- cbind(prAB_values_ly, prAB_value_ly)
        prAB_values_mle <- cbind(prAB_values_mle, prAB_value_mle)
      }
      
    }
    
  }
  
  if(loop == 1) {
    prAB_result <- prAB_values
    prAB_result_db <- prAB_values_db
    prAB_result_ly <- prAB_values_ly
    prAB_result_mle <- prAB_values_mle
    
    colnames(prAB_result) <- c("12", "13", "23")
    colnames(prAB_result_db) <- c("12", "13", "23")
    colnames(prAB_result_ly) <- c("12", "13", "23")
    colnames(prAB_result_mle) <- c("12", "13", "23")
    
  } else {
    prAB_result <- rbind(prAB_result, prAB_values)
    prAB_result_db <- rbind(prAB_result_db, prAB_values_db)
    prAB_result_ly <- rbind(prAB_result_ly, prAB_values_ly)
    prAB_result_mle <- rbind(prAB_result_mle, prAB_values_mle)
  }
  
  write.csv(prAB_result, 
            file=sprintf("%s/prAB_result.csv", folder_), 
            append=FALSE, 
            col.names=FALSE, 
            row.names=FALSE)
  
  write.csv(prAB_result_db, 
            file=sprintf("%s/prAB_result_dabrowska.csv", folder_), 
            append=FALSE, 
            col.names=FALSE, 
            row.names=FALSE)
  
  write.csv(prAB_result_ly, 
            file=sprintf("%s/prAB_result_linying.csv", folder_), 
            append=FALSE, 
            col.names=FALSE, 
            row.names=FALSE)
  
  write.csv(prAB_result_mle, 
            file=sprintf("%s/prAB_result_mle.csv", folder_), 
            append=FALSE, 
            col.names=FALSE, 
            row.names=FALSE)
  
}

