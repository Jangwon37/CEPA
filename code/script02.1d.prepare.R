# library
library(survival);
source('cOPT/R/opt1d.R');
source('cOPT/R/opt2d.R');
source('script99.misc.R');

opt1d.plot.recov = function(S,X,xlab='Time',ylab='Probability',main='') {
    par(cex.main=1.8)
	x = S[,1:2]; 
	r = range(x[is.finite(x)]);
	plot(survfit(X~1),conf.int=FALSE,xlim=c(r[1],53),ylim=c(0,1),xlab=xlab,ylab=ylab,main=main,cex.lab=1.5, lwd=2);
	x = 0; y = 1;
        S = rbind(c(0,S[1,1],0,0),S);
        for( i in 1:nrow(S) ) {
                x = c(x,S[i,2]);
                y = c(y,y[i]-S[i,3]*(S[i,2]-x[i]));
        }
	lines(x,y,type='l',col=4,lwd=2);
    
	x = S[,1:2]; 
	r = range(x[is.finite(x)]);
	x[is.infinite(x)] = r[2];
	y = cbind(S[,4],S[,4]);
	lines(as.vector(t(x)),as.vector(t(y)),type='l',col=3,lty=1,lwd=2);

	x = unique(as.vector(S[,1:2])); x = x[is.finite(x)];
	for( i in 1:length(x) ) {
		lines(c(x[i],x[i]),c(-100,100),col=2,lty=2);
	}
}

# read data
RECOV.TIME = read.csv('../data/censored/floor_disease.time_sep.csv',header=TRUE);
RECOV.EVENT = read.csv('../data/censored/disease.event_sep.csv',header=TRUE);

X = NULL;
for( i in 1:7 ) {
	x = Surv(RECOV.TIME[,i+1],RECOV.EVENT[,i+1]);
	if( i == 1 ) { X = x; }
	else { X = cbind(as.data.frame(X),x); }
}

colnames(X) = colnames(RECOV.TIME)[2:8];
row.names(X) = RECOV.TIME[,1];
CNX = colnames(X);
CNX.FULL = c('Infection','Eye','Respiratory','Digestive','Skin','Injury', 'Ear');

# opt 1d
S1D.ALL = vector('list',7);
for( i in 1:ncol(X) ) S1D.ALL[[i]] = opt1d(X[,i],A=c(1,53),max.depth=12,max.n=1,int.part=TRUE);
# plot
setEPS();
# postscript(file='fig.1d_surv_all.eps',width=9,height=9);
tiff(file='../result/censored/fig.1d_surv_all.tiff', width = 3000, height = 3000, units = "px", res=300);
par(mfrow=c(3,3),cex.axis=1.8, mar=c(4.5,4.5,3,1));
for( i in 1:ncol(X) ) opt1d.plot.recov(S1D.ALL[[i]],X[,i],main=CNX.FULL[i],xlab='Weeks After Diagnosis');
dev.off();


# median times
MED.TIMES = NULL;
for( i in 1:ncol(X) ) {
	s = S1D.ALL[[i]];
	t1 = icdf(s[,2],s[,4],0.5,given.occurrence=FALSE);
	t2 = icdf(s[,2],s[,4],0.5,given.occurrence=TRUE);
	#s = S1Dc.ALL[[i]];
	#t3 = icdf(s[,2],s[,4],0.5,given.occurrence=FALSE);
	#t4 = icdf(s[,2],s[,4],0.5,given.occurrence=TRUE);
	#MED.TIMES = rbind(MED.TIMES,c(t1,t2,t3,t4));
    MED.TIMES = rbind(MED.TIMES,c(t1,t2));
}
d.out = cbind(as.data.frame(CNX),round(10*MED.TIMES)/10);
#colnames(d.out) = c('event','med.all','med.obs','med.all.compl','med.obs.compl');
colnames(d.out) = c('event','med.all','med.obs');
write.table(d.out,file='../result/censored/tbl.median_time.txt',row.names=FALSE,quote=FALSE,sep='\t');

