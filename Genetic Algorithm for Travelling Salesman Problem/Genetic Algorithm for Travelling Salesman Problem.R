#生成随机地图
getMap = function(n,xlim=c(0,1),ylim=c(0,1))
{
	x = runif(n,xlim[1],xlim[2])
	y = runif(n,ylim[1],ylim[2])
	ans = as.matrix(dist(cbind(x,y)))
	return(list(x=x,y=y,ans=ans))
}

#生成单位圆地图
getRoundMap = function(n)
{
	theta = seq(0,2*pi,length.out=n+1)[1:n]
	x = cos(theta)
	y = sin(theta)
	ans = as.matrix(dist(cbind(x,y)))
	return(list(x=x,y=y,ans=ans))
}

#粗糙地计算总路程上界
UpperBound = function(disM)
{
	mx = apply(disM,1,max)
	return( sum(mx) )
}

#生成随机染色体
rndDNA = function(n)
{
	return( sample(1:n,n) )
}

#对一个解计算总路程的距离
calcScores = function(dna,disM)
{
	n = length(dna)
	tmp = cbind(dna[1:n],c(dna[2:n],dna[1]))
	len = apply(tmp,1,function(x) disM[x[1],x[2]])
	len = sum(len)
	return(len)
}

#根据每条染色体的分数计算权重，并以此抽样
roller = function(scores,k)
{
	scores = max(scores)-scores+1
	props = scores/sum(scores)
	N = length(scores)
	mxind = which.max(scores)#保留最优染色体
	ans = sample((1:N)[-mxind],k-1,replace = F,prob = props[-mxind])
	return(c(mxind,ans))
}

#种群中的繁殖过程
crossEvolve = function(i,nGroup,crossGroup,prop)
{
	a = nGroup[i,]
	b = crossGroup[i,]
	n = length(a)
	m = max(1,trunc(n*prop))

	tmpa = a

	st = sample(1:n,1)
	ind = st:(st+m)
	if (st+m>n)
	{
		bind = which(ind>n)
		ind[bind] = ind[bind] %% n + 1
	}
	cross = intersect(b,a[ind])
	tmpa[ind] = cross
	return(tmpa)
}

#染色体的自我变异
selfVariation = function(dna,prop)
{
	n = length(dna)
	pos = which(runif(n)<prop)
	if (length(pos)==0)
		return(dna)
	pos = sample(pos,1)
	newind = sample((1:n)[-pos],1)

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

#以某条染色体代表的解做图
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
	lines(x[dna[c(n,1)]],y[dna[c(n,1)]])
}

#结尾部分的动画函数，其中的具体参数需要变动
saveAnimation = function(species,dots,sleep=0.1)
{
	require(animation)
	scores = apply(species[[1]],1,calcScores,dots[[3]])
	n = length(scores)
	oopt = ani.options(interval = sleep, nmax = 10000, ani.width = 560, ani.height = 600)
	saveGIF(
	{
		for (i in 1:5)
			drawIt(dots,species[[1]][1,],xlab='',ylab='',main=scores[1],sub=paste('Generation:',species[[2]][1]))
		for (i in 1:9)
			drawIt(dots,species[[1]][i*20,],xlab='',ylab='',main=scores[i*20],sub=paste('Generation:',species[[2]][i*20]))
		for (i in 19:28)
				drawIt(dots,species[[1]][i*10,],xlab='',ylab='',main=scores[i*10],sub=paste('Generation:',species[[2]][i*10]))
		for (i in 290:410)
			if (i %% 3 == 0 )
				drawIt(dots,species[[1]][i,],xlab='',ylab='',main=scores[i],sub=paste('Generation:',species[[2]][i]))
		for (i in 411:443)
			if (i %% 2 == 0 )
				for (j in 1:2)
					drawIt(dots,species[[1]][i,],xlab='',ylab='',main=scores[i],sub=paste('Generation:',species[[2]][i]))
		for (i in 444:453)
			for (j in 1:5)
				drawIt(dots,species[[1]][i,],xlab='',ylab='',main=scores[i],sub=paste('Generation:',species[[2]][i]))
		for (i in 1:20)
			drawIt(dots,species[[1]][n,],xlab='',ylab='',main=scores[n],sub=paste('Generation:',species[[2]][n]))
	})
}

#Genetic Algorithm for Traveller Salesman Problem
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
	elite = mat.or.vec(maxElite,n)
	elitecount = 0
	eracount = 0

	LastEra = maxL
	outputRecorder = NULL
	GenerationRecorder = NULL
	showScore=FALSE
	eliteInto=FALSE
    #初始化结束
    
	while (!stopFlag)
	{
		cat('Generation:',iterCount,'Era:',eracount,'Elite:',elitecount)
		scores = apply(Group,1,calcScores,disM)
		bestScore = min(scores)
		mind = which.min(scores)
		bestDNA = Group[mind,]#记录最佳染色体

        #更新时代、精英的信息
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

		if (stayCount == maxStay)
		{
			stayCount = 0
			eliteBest = maxL
			elitecount = elitecount+1
			elite[elitecount,] = eliteDNA
			scores = apply(Group,1,calcScores,disM)
			a = which(scores == min(scores))
			nind = sample((1:N)[-a],length(a))
			Group[a,] = Group[nind,]
			eliteInto =TRUE
			scores = apply(Group,1,calcScores,disM)
			bestScore = min(scores)
			mind = which.min(scores)
			bestDNA = Group[mind,]
		}

		if (elitecount==maxElite)
		{
			Group[1:elitecount,]=elite
			elite = mat.or.vec(maxElite,n)
			elitecount = 0
			stayCount = 0
			eliteBest = maxL
			eracount = eracount+1
			maxStay = maxStay + 50
			showScore=TRUE
			scores = apply(Group,1,calcScores,disM)
			bestScore = min(scores)
			mind = which.min(scores)
			bestDNA = Group[mind,]
		}
        
        #对种群计算分数，产生繁殖与变异
		succind = roller(scores,N/2)
		nGroup = Group[succind,]
		crossind = sample(succind,N/2)
		crossGroup = Group[crossind,]
		crossans = t(sapply(1:(N/2),crossEvolve,nGroup,crossGroup,cp))
		crossGroup = rbind(nGroup,crossans)

		Group = t(apply(crossGroup,1,selfVariation,vp))
		if (eliteInto)
			eliteInto = FALSE
		else
			Group[1,] = bestDNA	

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
	if (drawing)
		drawIt(dots,allDNA)
	return(list(DNA = outputRecorder,Generation=GenerationRecorder))
}

roundots = getRoundMap(100)
system.time(( species = GA4TSP(dots=roundots,initDNA=NULL,N=50,cp=0.1,vp=0.01,maxIter=50000,maxStay=50,maxElite=2,drawing=TRUE) ))

#require(TSP)
#data("USCA50")
#tsp <- USCA50
#system.time(( species = GA4TSP(dots=list(a=NULL,b=NULL,disM=as.matrix(tsp)),initDNA=NULL,N=50,cp=0.1,vp=0.01,maxIter=50000,maxStay=50,maxElite=2,drawing=FALSE) ))
