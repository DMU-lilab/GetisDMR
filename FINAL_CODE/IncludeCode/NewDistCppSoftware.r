suppressMessages(library(geoR))
suppressMessages(library(inline))
suppressMessages(library(Rcpp))
suppressMessages(library(RcppArmadillo))
	
#Sys.setenv("CLINK_CPPFLAGS"="-I /usr/lib/R/site-library/RcppArmadillo/include/");
Sys.setenv("CLINK_CPPFLAGS"=ENVRCPPARMADILLO);
Sys.setenv("PKG_LIBS"="-llapack -lblas")
sourceCpp(MydistSoftware)
GetisNew<-function(chr,Z,OutCovFile)
{
	if(is.null(dim(chr))) chr=matrix(chr,ncol=1);
        chrsave=chr;
        Zsave=Z;
        if(nrow(chr)>10000)
        {
                chr=chr[1:10000,]
                Z=Z[1:10000]
        }

	if(is.null(dim(chr))) chr=matrix(chr,ncol=1);

        re=cbind(chr[,1],Z);
        breaks=seq(from=0,to=50000,by=50)
        v2=variog(coords=cbind(1:nrow(re),re[,1]),data=re[,2],breaks=breaks)
        start=500;
        RR=NULL;
        run=TRUE;
        num=0;
        step=50;
        while(run & num<max(v2$u))
        {
                x1=v2$u[v2$u<start+num]
                y1=v2$v[v2$u<start+num]
                RR=rbind(RR,c(start+num,summary(lm(y1~x1))$adj.r.squared))
                num=num+step;
#               if(RR[nrow(RR),2]<0.3) run=FALSE
        }

        a1=RR[1:(nrow(RR)-2),2]
        a2=RR[2:(nrow(RR)-1),2]
        a3=RR[3:(nrow(RR)-0),2]

        fs=apply(cbind(a1,a2,a3),1,mean)
        se=which(fs==max(fs))[1];
	if(se>80) se=80;
        as=c(RR[se,2],RR[se+1,2],RR[se+2,2])
#	print(as);
        add=(which(as==max(as))-1)[1];
#	print(add);
        se=se+add
#	print(se);

	
      	x1=v2$u[v2$u<RR[se,1]]
        y1=v2$v[v2$u<RR[se,1]]

        model=lm(y1~x1);
        vario=predict(model,data.frame("x1"=1:RR[se,1]))
        cov=(var(Z)-vario/2)/var(Z);
        pos=(1:RR[se,1])
        exc=(which(cov<0))
        if(length(exc)>0)
        {
                exc=min(exc)
                cov=cov[1:(exc-1)]
                pos=(1:RR[se,1])[1:(exc-1)]
        }
        pos=c(0,pos)
        cov=c(1,cov);

	write.table(cbind(pos,cov),OutCovFile,col.names=F,row.names=F,sep='\t');
        Z=Zsave;
        chr=chrsave;
        X_bar=mean(Z);
        S=sqrt(sum(Z^2)/length(Z)-X_bar^2)
        G.i=rep(NA,length(Z))

	Gnew=MyDist(Pos=matrix(chr[,1],ncol=1), Z=matrix(Z,ncol=1), COV=matrix(cov,ncol=1),X_bar, S, G=matrix(rep(0,length(G.i)),ncol=1),maxpos=max(pos))

	Gnew

}

