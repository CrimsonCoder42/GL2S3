#Matrix and vector calculations in Z5
#These functions use -1 and -2 aas names instead of 4 and 3
source("permutecalc.R")
Z5Sum <- function(a,b){
  s <-(a+b)%%5
  if (s > 2) s <- s-5
  return (s)
}

Z3Sum <- function(a,b) {
  s <- (a+b)%%3
  if (s > 1) s <- s-3
  return(s)
}

Z5Prod <- function(a,b){
  s <-(a*b)%%5
  if (s > 2) {
    s <- s-5
    }
  return (s)
}

Z3Prod <- function(a,b) {
  s <-(a*b)%%3
  if (s > 1) {
    s <- s-3
    }
  return(s)
}
Z3Prod(-1,0)


Z5Inv <- function(a){
  validate(need(a != 0, "Cannot divide by zero"))
  if (a == 1) return(1)
  if (a == 2) return(-2)
  if (a == -2) return(2)
  if (a == -1) return(-1)
}

Z3Inv <- function(a){
  validate(need(a != 0, "Cannot divide by zero"))
  if (a == 1) return(1)
  if (a == -1) return(-1)
}

#These vectorized functions make it possible 
#to multiply a matrix or a vector by a constant
vZ5Sum <- Vectorize(Z5Sum,c("a","b"));vZ5Sum 
vZ5Prod <- Vectorize(Z5Prod,c("a","b"));vZ5Prod

# Z3 Vectorize Functions
vZ3Sum <- Vectorize(Z3Sum, c("a","b"));vZ3Sum
vZ3Prod <- Vectorize(Z3Prod, c("a","b"));vZ3Prod

#5 rows for each of the 6 subspaces -- permits easy lookup
makeVecList <- function() {
  m <- matrix(ncol = 2, nrow = 30)
  v <- c(1,0)
  m[1,] <- vZ5Prod(0,v)
  m[2,] <- vZ5Prod(1,v)
  m[3,] <- vZ5Prod(2,v)
  m[4,] <- vZ5Prod(-2,v)
  m[5,] <- vZ5Prod(-1,v)
  v <- c(0,1)
  m[6,] <- vZ5Prod(0,v)
  m[7,] <- vZ5Prod(1,v)
  m[8,] <- vZ5Prod(2,v)
  m[9,] <- vZ5Prod(-2,v)
  m[10,] <- vZ5Prod(-1,v)
  v <- c(1,1)
  m[11,] <- vZ5Prod(0,v)
  m[12,] <- vZ5Prod(1,v)
  m[13,] <- vZ5Prod(2,v)
  m[14,] <- vZ5Prod(-2,v)
  m[15,] <- vZ5Prod(-1,v)
  v <- c(1,2)
  m[16,] <- vZ5Prod(0,v)
  m[17,] <- vZ5Prod(1,v)
  m[18,] <- vZ5Prod(2,v)
  m[19,] <- vZ5Prod(-2,v)
  m[20,] <- vZ5Prod(-1,v)
  v <- c(1,-2)
  m[21,] <- vZ5Prod(0,v)
  m[22,] <- vZ5Prod(1,v)
  m[23,] <- vZ5Prod(2,v)
  m[24,] <- vZ5Prod(-2,v)
  m[25,] <- vZ5Prod(-1,v)
  v <- c(1,-1)
  m[26,] <- vZ5Prod(0,v)
  m[27,] <- vZ5Prod(1,v)
  m[28,] <- vZ5Prod(2,v)
  m[29,] <- vZ5Prod(-2,v)
  m[30,] <- vZ5Prod(-1,v)
  return (m)
}
makeVecList()

#3 rows for each of the 4 subspaces -- permits easy lookup
makeVecListinZ3 <- function() {
  m <- matrix(ncol = 2, nrow = 12)
  v1 <- c(1,0)
  m[1,] <- vZ3Sum(0,v1)
  m[2,] <- vZ3Sum(1,v1)
  m[3,] <- vZ3Sum(-1,v1)
  v2 <- c(0,1)
  m[4,] <- vZ3Sum(0,v2)
  m[5,] <- vZ3Sum(1,v2)
  m[6,] <- vZ3Sum(-1,v2)
  v3 <- c(1,1)
  m[7,] <- vZ3Sum(0,v3)
  m[8,] <- vZ3Sum(1,v3)
  m[9,] <- vZ3Sum(-1,v3)
  v4 <- c(1,-1)
  m[10,] <- vZ3Sum(0,v4)
  m[11,] <- vZ3Sum(1,v4)
  m[12,] <- vZ3Sum(-1,v4)
  return (m)
}
makeVecListinZ3()

library(prodlim) #needed for row.match
#Multiply a matrix by a vector in Z5
ActOnVector <- function(A,v){
  return(vZ5Sum(vZ5Prod(v[1],A[,1]),vZ5Prod(v[2],A[,2])))
}

ActOnVectorinZ3 <- function(A,v){
  return(vZ3Sum(vZ3Prod(v[1],A[,1]),vZ3Prod(v[2],A[,2])))
}

#Applies matrix A to a vector from subspace idx
#Looks up the result to find what subspace it is in.
Transform <- function(A,idx){
  m <- makeVecList()
  v <- m[5*idx-3,]   #second vector (first nonzero one) in subspace idx
  x <- vZ5Sum(vZ5Prod(v[1],A[,1]),vZ5Prod(v[2],A[,2]))
#Find what row in matrix m amatches vector v
  r <- row.match(x,m)
  return(floor((r+4)/5))
}

TransforminZ3 <- function(A,idx){
  m <- makeVecListinZ3()
  v <- m[2*idx-3,]   #second vector (first nonzero one) in subspace idx
  x <- vZ3Sum(vZ3Prod(v[1],A[,1]),vZ3Prod(v[2],A[,2]))
  #Find what row in matrix m amatches vector v
  r <- row.match(x,m)
  return(floor((r+4)/5))
}

# let <- TransforminZ3(testMatrix,6);let
# 
# testMatrix <- matrix(c(0,1,1,0),2); testMatrix
# test <- TransforminZ3(testMatrix, 1);test
# 
# fval <- sapply(1:6,TransforminZ3,A=testMatrix); fval
# test <- Perm.cycle.convert(fval);test

#Multiply 2x2 matrices in Z5
Z3MatProd <- function(A,B) {
  v1 <- B[,1]
  x1 <- vZ3Sum(vZ3Prod(v1[1],A[,1]),vZ3Prod(v1[2],A[,2]))
  v2 <- B[,2]
  x2 <- vZ3Sum(vZ3Prod(v2[1],A[,1]),vZ3Prod(v2[2],A[,2]))
  return(cbind(x1,x2))
}

#Multiply 2x2 matrices in Z3
Z5MatProd <- function(A,B) {
  v1 <- B[,1]
  x1 <- vZ5Sum(vZ5Prod(v1[1],A[,1]),vZ5Prod(v1[2],A[,2]))
  v2 <- B[,2]
  x2 <- vZ5Sum(vZ5Prod(v2[1],A[,1]),vZ5Prod(v2[2],A[,2]))
  return(cbind(x1,x2))
}

Z5Det <- function(A) {
  Z5Sum(Z5Prod(A[1,1],A[2,2]),Z5Prod(A[1,2],A[2,1]))
}

#Create a random matrix with specified determinant and trace
Z5CreateMatrix <- function(det,trc){
  elements <- c(0,1,2,-2,-1)
#generate random top left entry
  a11 <- as.numeric(sample(elements,1))
#make the trace correct
  a22 <- Z5Sum(-a11,trc)
#choose a random nonzero element
  b1 <- sample(elements[2:5],1)   #nonzero
  #Required product for the off-diagonal elements
  offDiag <- Z5Sum(-det,Z5Prod(a11,a22))    #could be zero
  #Multiply it by the inverse of the known element
  b2 <- Z5Prod(offDiag,Z5Inv(b1))     #could be zero
  #The following makes it possible to have a zero in either off-diagonal position
  if (sample(2,1)==1)
    return (matrix(c(a11,b1,b2,a22),2))
  return (matrix(c(a11,b2,b1,a22),2))
}

Z3CreateMatrix <- function(det,trc){
  elements <- c(0,1,-1)
  #generate random top left entry
  a11 <- as.numeric(sample(elements,1))
  #make the trace correct
  a22 <- Z3Sum(-a11,trc)
  #choose a random nonzero element
  b1 <- sample(c(-1,1),1)   #nonzero
  #Required product for the off-diagonal elements
  offDiag <- Z3Sum(-det,Z3Prod(a11,a22))    #could be zero
  #Multiply it by the inverse of the known element
  b2 <- Z3Prod(offDiag,Z5Inv(b1))     #could be zero
  #The following makes it possible to have a zero in either off-diagonal position
  if (sample(2,1)==1)
    return (matrix(c(a11,b1,b2,a22),2))
  return (matrix(c(a11,b2,b1,a22),2))
}

A <- Z3CreateMatrix(1,0);A
