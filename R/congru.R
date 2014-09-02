congru <- function(x,y){
  x=c(x); y=c(y)
  if(length(x)!=length(y)){stop("x and y must have same length")}
  sum(x*y)/sqrt(sum(x^2)*sum(y^2))
}