#!/usr/bin/env Rscript
# Rscript to get the elbow point of a curve
# Ref: http://stackoverflow.com/questions/24895916/vector-projection-to-determine-elbow-in-fitted-curve
args <- commandArgs(trailingOnly = TRUE)
input=args[1]

get_elbow_index <- function(x, y){
  
  #fake data
  #x <- 5:300
  #y <- (x - 0.03*x^2 + 0.02*x^3 + rnorm(length(x), sd=5000))/1000
  
  myData <- data.frame(x, y)
  
  # fitted curve (I used a simpler example)
  result <- lm(y ~ x + I(x^2) + I(x^3) , data=myData)
  p <- fitted(result)
  
  # line connecting endpoints of fitted curve
  i1 <- which.min(x)
  i2 <- which.max(x)
  slope <- (p[i2] - p[i1]) / (x[i2] - x[i1])
  int <- p[i1] - slope*x[i1]
  
  # for every point on the predicted curve (xi, pi), the perpendicular line that goes through that point has
  perpslope <- -1/slope
  perpint <- p - perpslope*x
  
  # the intersection of the perp line(s) with the connecting line is
  xcross <- (int - perpint) / (perpslope - slope)
  ycross <- slope*xcross + int
  
  # the distance between the intersection and the point(s) is
  dists <- sqrt((x - xcross)^2 + (y - ycross)^2)
  
  # the index of the farthest point
  elbowi <- which.max(as.vector(dists))
  
  return(elbowi)
  
}


df=read.table(input, header=F)$V1
Fn=ecdf(log10(df))
pdf(gsub(".txt",".pdf",input), width=5, height=5)
plot(Fn, xaxt="n", ylab="Cummulative percentage", xlab="Distance between neighboring HTNEs (log10)")
y1=floor(range(knots(Fn)))
pow <- seq(y1[1], y1[2]+1)
ticksat <- as.vector(sapply(pow, function(p) log10((1:9)*10^p)))
axis(1, pow, labels = 10^pow, las=1)
axis(1, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
x=seq(range(knots(Fn))[1], range(knots(Fn))[2],by=(range(knots(Fn))[2]-range(knots(Fn))[1])/500); y=Fn(x);
elbowi = get_elbow_index(x,y) # the point with the farthest distance
xi=x[elbowi]; # == 4.281637
inv_ecdf <- function(f){ x <- environment(f)$x; y <- environment(f)$y; approxfun(y, x)}; g <- inv_ecdf(Fn);
xi=g(0.9)
#xi=log10(50000)
lines(x=c(1,xi,xi),y=c(Fn(xi),Fn(xi),0), col='red'); 
points(x=xi,y=Fn(xi),col='red'); 
text(x=xi,y=Fn(xi),labels=paste0(" ", round(10**xi,0), ", ", round(100*Fn(xi),2),"%"), adj=c(0,1), col='red')
dev.off()



# # plot the data
# eqscplot(x, y)
# lines(x[c(i1, i2)], p[c(i1, i2)])
# points(x[elbowi], p[elbowi], pch=16, col="red")
# lines(x[order(x)], p[order(x)], col="blue")
# lines(c(x[elbowi], xcross[elbowi]), c(p[elbowi], ycross[elbowi]), col="red")