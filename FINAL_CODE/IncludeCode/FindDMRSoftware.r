# Fix Small Mistakes, adding the case where the first or the last are FALSE, but in the middle it meets the criteria #
# Copied from regionAsBedFeb22
## Further fix mistakes from the original codes ###
# Fix Small Mistakes, adding the case where the first or the last are FALSE, but in the middle it meets the criteria #
regionAsBed <- function(marker, cf.length=5, tolerance.length=NULL,cut=NULL,chrom)
{
	r <- rle(marker);
	markersave=marker;
	
	if(is.null(tolerance.length))
	{
	        ##The tolerance length is NUL, which means any number of false
	        ##values, even 1, can break a regionn into two.
	        end <- with(r,cumsum(lengths)[values & lengths>cf.length])
	        start <- end - with(r,lengths[values & lengths>cf.length]) + 1
	        #df <- data.frame(chrom,start,end)
	        if(length(start) == 0 || length(end) == 0)
		{
            		df <- data.frame(chrom=character(), start=integer(), end=integer())
	        } 
		else 
		{
			df <- data.frame(chrom, start, end)
        	}
	} 
	else 
	{
		# These are the situations where the start or the end are FALSE #
   		end1 <- with(r,cumsum(lengths)[values & lengths>cf.length-tolerance.length])
	        start1 <- end1 - with(r,lengths[values & lengths>cf.length-tolerance.length]) + 1
	        ##Tolerance length is not null, which means if the number of false
	        ##is less than or equal to tolerance length, the two regions will
	        ##be combined together
	        ##Get all the starts and ends
	        end<-cumsum(r$lengths)
	        start<-end - r$lengths + 1
	        revision.mk <- with(r, lengths<=tolerance.length & !values)
		##The position that need to be revised
	        revision.posi <- data.frame(cbind(start,end))[revision.mk,]
	        ##revise the orignal marker 
	        invisible(apply(revision.posi, 1, function(x) marker[x[1]:x[2]] <<- T))
	        r <- rle(marker)
		end <- with(r,cumsum(lengths)[values & lengths>cf.length])
	        start <- end - with(r,lengths[values & lengths>cf.length]) + 1;
		
		if(length(start)>0)
		{
			
			inc=rep(TRUE,length(start))
			for(i in 1:length(inc))
			{
				tmp=mean(markersave[start[i]:end[i]]);
				tmp1=(1-tolerance.length/cf.length)
				inc[i]= tmp>tmp1;
			}
			start=start[inc];
			end=end[inc];
		}
		
		start=c(start,start1);
		end=c(end,end1);
		
	        if(length(start) == 0 || length(end) == 0)
		{
			df <- data.frame(chrom=character(), start=integer(), end=integer())
		} 
		else 
		{
			df <- data.frame(chrom, start, end)
		}
	}
	return(df)
}




findDMR<-function(G.i,posi,range.di.high,range.di.low=NULL,cf.length=5, tolerance.length=NULL,outfile,chrom)
{
	if(is.null(range.di.low)) range.di.low=-range.di.high;
	if(length(range.di.low)!=length(range.di.high)) stop("The lenght of lower bound does not equal to the high bound, please check");

	########## FOR CURRENT USE ONLY, WILL MODIFY ONCE A DECISION HAS BEEN MADE #############

	for(j in 1:length(G.i))
	{
		g.i=G.i[[j]];
		re=list();
		re$Gi=g.i;
		for(jj in 1:length(range.di.high))
		{
			#cat(range.di[jj]);
			cf.di=c(range.di.high[jj],range.di.low[jj])
			dffi1=round(cf.di[1],2)
			re$mk <- re$Gi > cf.di[1];

			re$posi=posi;
			out=list();
			flag<-0;
			bigout=TRUE;
			bd<-tryCatch(regionAsBed(marker=as.vector(re$mk), cf.length=cf.length, tolerance.length=tolerance.length, chrom=chrom),error = function(e) {flag<<-1})
			if(flag == 1 || nrow(bd) == 0) 
			{
      				#message("error/No big in processing large deletion for ",chrom)
				bigout=FALSE;	
			} 		
			else 
			{
				bd.onchrom<-apply(bd, 1, function(x){s<-as.integer(x[2]);e<-as.integer(x[3]);c(re$posi[s],re$posi[e])})
				bd.onchrom <- t(bd.onchrom)
				bd.onchrom <- data.frame(chrom, bd.onchrom)
				colnames(bd.onchrom) <- c("chrom","start","end")
				out$large <- bd.onchrom
			}

			dffi2=round(cf.di[2],2)
			re$mk <- re$Gi < (cf.di[2])
			flag<-0
			smallout=TRUE;
			bd<-tryCatch(regionAsBed(marker=as.vector(re$mk), cf.length=cf.length, tolerance.length=tolerance.length, chrom=chrom),error = function(e) {flag<<-1})
			if(flag == 1 || nrow(bd) == 0) 
			{
				#message("error/No small in processing small deletion for ",chrom)
				smallout=FALSE;
			} 
			else 
			{
				bd.onchrom<-apply(bd, 1, function(x){s<-as.integer(x[2]);e<-as.integer(x[3]);c(re$posi[s],re$posi[e])})
				bd.onchrom <- t(bd.onchrom)
				bd.onchrom <- data.frame(chrom, bd.onchrom)
				colnames(bd.onchrom) <- c("chrom","start","end")
				out$small <- bd.onchrom
			}

			if(bigout) write.table(with(out$large, paste(chrom, ":", start, ":", end, "   ", end - start, sep="")),file=paste(outfile,"G",j,"D",dffi1,sep=""), row.name=F, col.name=F, quote=F, sep="\t")
			if(smallout) write.table(with(out$small, paste(chrom, ":", start, ":", end, "   ", end - start, sep="")),file=paste(outfile,"G",j,"D",dffi2,sep=""), row.name=F, col.name=F, quote=F, sep="\t")
		}
	}
}





