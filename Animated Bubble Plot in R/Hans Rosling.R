require(RJSONIO)
dat = fromJSON('Wealth and Health.txt')

n = length(dat)
tbls2 = NULL

for (i in 1:n) {
    dti = dat[[i]]
    dfr = initdfr
    dfr$name = rep(dti$name, 210)
    dfr$region = rep(dti$region, 210)
    for (inc in dti$income) {
        ii = which(dfr$year == as.integer(inc[1]))
        dfr$income[ii] = inc[2]
    }
    for (pop in dti$population) {
        ii = which(dfr$year == as.integer(pop[1]))
        dfr$pop[ii] = pop[2]
    }
    for (life in dti$lifeExpectancy) {
        ii = which(dfr$year == as.integer(life[1]))
        dfr$life[ii] = life[2]
    }

    tbls2 = rbind(tbls2, dfr)
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
