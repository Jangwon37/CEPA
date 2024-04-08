# data on the uniform distribution
random_data_unif = function(N,min,max,of=0) {
  # X: multivariate survival times
  # T: true event times
  if(length(min)!=length(max)) {print("Length must be same"); return();}
  X = NULL;  T = NULL; colxn = NULL; coltn = NULL;
  for(i in 1:length(min)){
    tmpT = runif(N,min[i],max[i]);
    tmpC = runif(N,min[1],max[length(max)])+of;
    tmpX = pmin(tmpT,tmpC); tmpD = (tmpT<=tmpC);
    X = cbind(X, tmpX, tmpD); T = cbind(T, tmpT);
    colxn = c(colxn, sprintf("X%d",i), sprintf("D%d",i));
    coltn = c(coltn, sprintf("T%d",i));
  }
  colnames(X) = colxn; colnames(T) = coltn;
  return( list(X=X, T=T) );
}

# data on the exponential distribution
random_data_exp = function(N,lt,lc,of=0) {
  # X: multivariate survival times
  # T: true event times
  X = NULL;  T = NULL; colxn = NULL; coltn = NULL;
  for(i in 1:length(lt)){
    tmpT = rexp(N,lt[i]);
    tmpC = rexp(N,lc[i])+of;
    tmpX = pmin(tmpT,tmpC); tmpD = (tmpT<=tmpC);
    X = cbind(X, tmpX, tmpD); T = cbind(T, tmpT);
    colxn = c(colxn, sprintf("X%d",i), sprintf("D%d",i));
    coltn = c(coltn, sprintf("T%d",i));
  }
  colnames(X) = colxn; colnames(T) = coltn;
  return( list(X=X, T=T) );
}

# clayton q=1
random_data_clayton = function(N,lambda=1/2,q=1,of=0) {
  U1 = runif(N); U2 = runif(N);
  a = (1-U2)^(-1/q);
  
  T1 = q*log( (1-a) + a*(1-U1)^(-1/(1+q)) );
  T2 = -log(1-U2);
  
  C1 = rexp(N,lambda)+of;
  C2 = rexp(N,lambda)+of;
  
  X1 = Surv(apply(cbind(T1,C1),1,min),(T1<=C1) );
  X2 = Surv(apply(cbind(T2,C2),1,min),(T2<=C2) );
  X = cbind(as.data.frame(X1),X2); colnames(X) = c('X1','X2');
  
  return( list(X=X, T=cbind(T1,T2)) );
}

PrAB = function(dist, xbound, ybound){
  tmp = rep(0,nrow(dist));
  
  # Pr[A<=B]
  P = 0;
  for(i in 1:nrow(dist)){
    x1 = dist[i,1]; x2 = dist[i,2];
    y1 = dist[i,3]; y2 = dist[i,4];
    p = dist[i,6];
    
    if(p<0||p>1||is.nan(p)){
      p=0;
    }
    if(x2<y1){
      P = P + p;
      tmp[i] = p;
    }
    else if(x2==Inf&&y2==Inf){
      P = P + p/2;
    }
    else if(y2==Inf){
      P = P + p;
    }
    else if(x2==Inf){
      
    }
    else if(x1<=y1&&y1<=x2&&x1<=y2&&y2<=x2){
      if((x1==x2&&y1==y2)||(x2==Inf)||(y2==Inf)){
        P = P + p;
        tmp[i] = p;
      }
      else{
        P = P + p*(((y1-x1)+(y2-x1))*(y2-y1)/2)/((x2-x1)*(y2-y1)); 
        tmp[i] = p*(((y1-x1)+(y2-x1))*(y2-y1)/2)/((x2-x1)*(y2-y1));
      }
    }
    else if(y1<=x1&&x1<=y2&&y1<=x2&&x2<=y2){
      if((x1==x2&&y1==y2)||(x2==Inf)||(y2==Inf)){
        P = P + p;
        tmp[i] = p;
      }
      else{
        P = P + p*(((y2-x1)+(y2-x2))*(x2-x1)/2)/((x2-x1)*(y2-y1)); 
        tmp[i] = p*(((y2-x1)+(y2-x2))*(x2-x1)/2)/((x2-x1)*(y2-y1))
      }
    }
    else if(y1<=x1&&x1<=y2&&x1<=y2&&y2<=x2){
      if((x1==x2&&y1==y2)||(x2==Inf)||(y2==Inf)){
        P = P + p;
        tmp[i] = p;
      }
      else{
        P = P + p*((y2-x1)*(y2-x1)/2)/((x2-x1)*(y2-y1));
        tmp[i] = p*((y2-x1)*(y2-x1)/2)/((x2-x1)*(y2-y1))
      }
    }
    else if(x1<=y1&&y1<=x2&&y1<=x2&&x2<=y2){
      if((x1==x2&&y1==y2)||(x2==Inf)||(y2==Inf)){
        P = P + p;
        tmp[i] = p;
      }
      else{
        P = P + p*(((x2-x1)*(y2-y1))-((x2-y1)*(x2-y1))/2)/((x2-x1)*(y2-y1));
        tmp[i] = p*(((x2-x1)*(y2-y1))-((x2-y1)*(x2-y1))/2)/((x2-x1)*(y2-y1))
      }
    }
  }
  
  return (list(P=P, sub.region = tmp));
}























