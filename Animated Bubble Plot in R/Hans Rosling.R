#Import the Strange Dataset
rawdata=readLines("Wealth and Health.txt")

#Get the useful words
data=unlist(strsplit(rawdata,"\\[\\{\"|\":\"|\",\"|\":\\[\\[|\\],\\[|\\]\\],\"|\\]\\]\\},\\{\"|\\]\\]\\}\\]"))
data=data[which(data!="")]
n=length(data)

#Set index
nameind=which(data=="name")
regind=which(data=="region")
incind=which(data=="income")
popind=which(data=="population")
lifind=which(data=="lifeExpectancy")
endpoint=nameind-1
endpoint=endpoint[2:180]
endpoint[180]=45999

name=data[nameind+1]
region=data[regind+1]

#An empty data.frame type
initdfr=data.frame(name=rep(0,210),region=rep(0,210),year=1800:2009,income=rep(0,210),pop=rep(0,210),life=rep(0,210))

#An string split function only working for the comma
splt=function(x) return(as.numeric(unlist(strsplit(x,","))))

#Linear interpolation function, for those zero data
itpl=function(a)
{
	ind=which(a>0)
	if (ind[1]>1)
		a[1:(ind[1]-1)]=rep(a[ind[1]],ind[1]-1)
	n=length(ind)
	if (ind[n]<length(a))
		a[(ind[n]+1):length(a)]=rep(a[ind[n]],length(a)-ind[n])
	for (i in 1:(n-1))
		a[ind[i]:ind[i+1]]=rep(a[ind[i]],ind[i+1]-ind[i]+1)+(seq(ind[i],ind[i+1],1)-ind[i])*(a[ind[i+1]]-a[ind[i]])/(ind[i+1]-ind[i])
	return(a)
}

#Drag data information from words
tbls=NULL
for (i in 1:180)
{
	dfr=initdfr
	dfr$name=rep(name[i],210)
	dfr$region=rep(region[i],210)
	for (j in (incind[i]+1):(popind[i]-1))
	{
		tmp=splt(data[j])
		ii=which(dfr$year==tmp[1])
		dfr$income[ii]=tmp[2]
	}
	for (j in (popind[i]+1):(lifind[i]-1))
	{
		tmp=splt(data[j])
		ii=which(dfr$year==tmp[1])
		dfr$pop[ii]=tmp[2]
	}
	for (j in (lifind[i]+1):endpoint[i])
	{
		tmp=splt(data[j])
		ii=which(dfr$year==tmp[1])
		dfr$life[ii]=tmp[2]
	}
	
	tbls=rbind(tbls,dfr)
}

#Two country with only one record, meaningless
ind=which(tbls$name=="Mayotte")
tbls=tbls[-ind,]
ind=which(tbls$name=="Tokelau")
tbls=tbls[-ind,]
name=name[c(-28,-177)]
region=region[c(-28,-177)]

#Linear interpolation
for (i in 1:178)
{
	ind=(210*(i-1)+1):(i*210)
	tbls$income[ind]=itpl(tbls$income[ind])
	tbls$pop[ind]=itpl(tbls$pop[ind])
	tbls$life[ind]=itpl(tbls$life[ind])
}

require(ggplot2)

#Draw function with ggplot2
drawit=function(yr,scl=15)
{
    ind=which(tbls$year==yr)
	d.f=data.frame(yr=yr)
    p=ggplot(aes(x=log(income),y=life,size=pop,colour=as.factor(region)),pch=21,data=tbls[ind,])
    p+geom_point(show_guide = FALSE)+
	geom_point(shape = 1,colour = "black",show_guide = FALSE)+
	xlim(5.5,11.7)+ylim(10,83)+scale_area(range = c(1, scl))+
	annotate("text", x=10, y=15, label = yr,size=30,color="grey")
}
#drawit(1800)

#Automatically repeat the drawing procedure
finaldraw=function(a,b)
{
	for (i in 1:10)
		print(drawit(a))
	for (i in a:b)
		print(drawit(i))
	for (i in 1:10)
		print(drawit(b))
}

#finaldraw(1800,2009)


require(animation)

#sett ffmpeg in Windows = =||
oopts = ani.options(ffmpeg = "D:/ffmpeg/bin/ffmpeg.exe")

#Use the function from animation to make the final movie
saveVideo({
    finaldraw(1800,2009)
	ani.options(interval = 0.1, nmax = 230)
}, video.name = "HansRosling.mp4", other.opts = "-b 500k")