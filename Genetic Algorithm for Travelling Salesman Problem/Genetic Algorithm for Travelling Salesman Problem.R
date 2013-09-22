getMap = function(n,xlim=c(0,1),ylim=c(0,1))
{
	x = runif(n,xlim[1],xlim[2])
	y = runif(n,ylim[1],ylim[2])
	ans = as.matrix(dist(cbind(x,y)))
	return(list(x=x,y=y,ans=ans))
}

getRoundMap = function(n)
{
	theta = seq(0,2*pi,length.out=n+1)[1:n]
	x = cos(theta)
	y = sin(theta)
	ans = as.matrix(dist(cbind(x,y)))
	return(list(x=x,y=y,ans=ans))
}

UpperBound = function(disM)
{
	mx = apply(disM,1,max)
	return( sum(mx) )
}

rndDNA = function(n)
{
	return( sample(1:n,n) )
}

calcScores = function(dna,disM)
{
	n = length(dna)
	tmp = cbind(dna[1:(n-1)],dna[2:n])
	len = apply(tmp,1,function(x) disM[x[1],x[2]])
	len = sum(len)
	return(len)
}

roller = function(scores,k)
{
	scores = max(scores)-scores+1
	props = scores/sum(scores)
	N = length(scores)
	mxind = which.max(scores)
	ans = sample((1:N)[-mxind],k-1,replace = F,prob = props[-mxind])
	return(c(mxind,ans))
}

crossEvolve = function(i,nGroup,crossGroup,prop)#random points
{
	a = nGroup[i,]
	b = crossGroup[i,]
	n = length(a)
	m = max(1,trunc(n*prop))

	tmpa = a

	st = sample(1:n,1)
	#ind = sample(1:n,m)
	ind = st:min(n,st+m)
	cross = intersect(b,a[ind])
	tmpa[ind] = cross

	return(tmpa)
}

selfVariation = function(dna,prop)
{
	n = length(dna)
	pos = which(runif(n)<prop)
	if (length(pos)==0)
		return(dna)
    pos = sample(pos,1)
	newind = sample((1:n)[-pos],1)
	#dna[c(pos,newind)] = dna[c(newind,pos)]
	if (pos>newind)
	{
		tmp = dna[newind:n]
		tmp = tmp[-(pos-newind+1)]
		dna[newind] = dna[pos]
		dna[(newind+1):n] = tmp
	}
	else
	{
		tmp = dna[1:newind]
		tmp = tmp[-pos]
		dna[newind] = dna[pos]
		dna[1:(newind-1)] = tmp
	}
	return(dna)
}

drawIt = function(dots,dna,xlab=NULL,ylab=NULL,main=NULL,sub=NULL)
{
	x = dots[[1]]
	y = dots[[2]]
	if (is.null(main))
	{
		scores = calcScores(dna,dots[[3]])
		plot(x,y,main=as.character(scores))
	}
	else
		plot(x,y,main=main,xlab=xlab,ylab=ylab,sub=sub)
	n = length(dna)
	for (i in 1:(n-1)) lines(x[dna[i:(i+1)]],y[dna[i:(i+1)]])
}

saveAnimation = function(species,dots,sleep)
{
	require(animation)
	scores = apply(species[[1]],1,calcScores,dots[[3]])
	n = length(scores)
	oopt = ani.options(interval = sleep, nmax = n+1, ani.width = 560, ani.height = 600)
	saveGIF(
	{
		for (i in 1:5)
			drawIt(dots,species[[1]][1,],xlab='',ylab='',main=scores[1],sub=paste('Generation:',species[[2]][1]))
		for (i in 1:n)
			drawIt(dots,species[[1]][i,],xlab='',ylab='',main=scores[i],sub=paste('Generation:',species[[2]][i]))
		for (i in 1:10)
			drawIt(dots,species[[1]][n,],xlab='',ylab='',main=scores[n],sub=paste('Generation:',species[[2]][n]))
	})
}

GA4TSP = function(dots,initDNA=NULL,N,cp,vp,maxIter,maxStay,maxElite,drawing)
{
	disM = dots[[3]]
	n = nrow(disM)
	if (N %%2 >0)
		N = N+1
	Group = t(sapply(rep(n,N),rndDNA))
	if (!is.null(initDNA))
	{
		if (is.null(initDNA))
			Group[1,]=initDNA
		else
		{
			if (nrow(initDNA)>N)
				initDNA = initDNA[1:N,]
			Group[1:nrow(initDNA),]=initDNA
		}
	}
	maxL = UpperBound(disM)
	stopFlag = FALSE
	iterCount = 1
	stayCount = 0
	allBest = maxL
	eliteBest = maxL
	elite = matrix(numeric(maxElite*n),nrow=maxElite)
	elitecount = 0
	eracount = 0

	LastEra = maxL
	outputRecorder = NULL
	GenerationRecorder = NULL
	showScore=FALSE

	while (!stopFlag)
	{
		cat('Generation:',iterCount,'Era:',eracount,'Elite:',elitecount)
		scores = apply(Group,1,calcScores,disM)
		succind = roller(scores,N/2)
		bestScore = min(scores)
		mind = which.min(scores)
		bestDNA = Group[mind,]
		#maxL = (max(scores)-allBest)*shrink+allBest

		nGroup = Group[succind,]
		crossind = sample(succind,N/2)
		crossGroup = Group[crossind,]
		crossans = t(sapply(1:(N/2),crossEvolve,nGroup,crossGroup,cp))
		crossGroup = rbind(nGroup,crossans)

		#Group = t(apply(crossGroup,1,selfVariation,vp))
		Group = t(apply(crossGroup,1,selfVariation,vp))
		Group[1,] = bestDNA

		if (bestScore<eliteBest)
		{
			stayCount = 0
			eliteBest = bestScore
			eliteDNA = bestDNA
			if (eliteBest<allBest)
			{
				allBest = eliteBest
				allDNA = eliteDNA
				outputRecorder = rbind(outputRecorder,allDNA)
				GenerationRecorder = c(GenerationRecorder,iterCount)
				if (drawing)
					drawIt(dots,allDNA,main=as.character(allBest))
			}
		}
		else
			stayCount = stayCount+1

		if (stayCount > maxStay)
		{
			stayCount = 0
			eliteBest = maxL
			elitecount = elitecount+1
			elite[elitecount,]=bestDNA
			scores = apply(Group,1,calcScores,disM)
			a = which(scores == max(scores))
			for (i in a)
				Group[a,] = rndDNA(n)
				#Group[which.max(scores)[1],] = rndDNA(n)
		}
		#if (elitecount==N/2)
		if (elitecount==maxElite)
		{
			Group[1:elitecount,]=elite
			elite = matrix(numeric(maxElite*n),nrow=maxElite)
			elitecount = 0
			stayCount = 0
			eliteBest = maxL
			eracount = eracount+1
			maxStay = maxStay + 50
			showScore=TRUE
		}

		stopFlag = (iterCount>=maxIter) 
		iterCount = iterCount+1
		cat(' Best:',bestScore,'All:',allBest,'Stay:',paste(stayCount,'/',maxStay,sep=''),'\n')
		if (showScore)
		{
			scores = apply(Group,1,calcScores,disM)
			show(scores[1:(N/2)])
			show(scores[(N/2+1):(N)])
			Sys.sleep(5)
			showScore=FALSE
		}
	}
	drawIt(dots,allDNA)
	return(list(DNA = outputRecorder,Generation=GenerationRecorder))
}

#dots = getMap(n=100,xlim=c(0,1),ylim=c(0,1))
#disM = dots[[3]]
#load('roundot.rda')
#dots = getMap(100)
roundots = getRoundMap(100)
system.time(( species = GA4TSP(dots=roundots,initDNA=NULL,N=50,cp=0.1,vp=0.01,maxIter=300000,maxStay=50,maxElite=2,drawing=FALSE) ))

save(species,file='species.rda')
