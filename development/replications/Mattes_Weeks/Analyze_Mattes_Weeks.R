## Table in paper printed at end of run!
rm(list=ls())

source('../Code/HOE_APSR.R')

options(device="quartz")
load('../Mattes_Weeks/MattesHOE.Rda')




p1<-pos_measure(output$h1,trim=0.025)
p2<-pos_measure(output$h2,trim=0.025)


plot(make.rank(p1$diff),p1$diff,type="l",ylim=range(c(p1$diff,p2$diff)),col="red",lwd=2)
lines(make.rank(p2$diff),p2$diff,type="l",lwd=3)
abline(h=0)

#lines(make.rank(p3$observed),p3$diff,type="l",col="blue")

gg_color <- function(n,alpha0=1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100,alpha=alpha0)[1:n]
}

pdf("../Mattes_Weeks/posmattesweeks.pdf",h=4*1.1,w=7*1.5/2)


cols<-gg_color(2)

par(mar = c(3, 3,1,0.2), # Dist' from plot to side of page
    mgp = c(2, 0.4, 0), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.01, # Reduce tick length
    xaxs = "i", yaxs = "i", # Remove plot padding
    oma =c(.2,.2,0,0))



plot(make.rank(p1$diff),p1$diff,type="n",col="red",lwd=2,ylim=range(c(p1$diff,p2$diff,0)),
     xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(-5,105))
#abline(c(0,1),col="gray50")
abline(h=0,lwd=2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray90")
abline(h=(-100:100)/2,col="white")
abline(v=(0:100)*10,col="white")

mtext(side=1,line=1.7,font=2,text="Percentile",las=0)
mtext(side=2,line=2,font=2,text="Excess Kurtosis",las=0)
axis(1)
axis(2)
lines(make.rank(p1$diff),p1$diff,col=cols[1],lwd=2)
lines(make.rank(p2$diff),p2$diff,col=cols[2],lwd=2)


legend.txt<-c("Hawks Experiment","Doves Experiment")
legend("topright",legend=legend.txt,col=gg_color(2),lty=1,xpd=T,lwd=2,bg="gray90",bty="n",box.col="gray90")

dev.off()


