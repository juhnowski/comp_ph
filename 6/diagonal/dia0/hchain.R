nn=8
nst=2^nn

mat<-matrix(0.0, nrow = nst, ncol = nst)
vec<-matrix(0.0, nrow = nst, ncol = nst)
enr<-matrix(0.0, nrow = 1, ncol = nst)
spn<-matrix(0.0, nrow = 1, ncol = nst)
mag<-matrix(0.0, nrow = 1, ncol = nst)

btest <- function (value, pos) {
  if (bitwAnd(value,bitwShiftL(1,pos))>0) {
    return(T)
  }
  return(F)
}

hamiltonian <- function(nn, nst, mat) {
# Hamiltonian
  result <- mat

  for (a in 1:nst-1) {
    for (i in 1:nn-1) {
      j=(i+1)%%nn
#      print(paste("a=",a ," , i=",i,", j=", j))
      if (btest(a,i)==btest(a,j)){
        mat[a+1,a+1]=mat[a+1,a+1]+0.25
#        print(paste("mat(a,a)=", mat[a+1,a+1], btest(a,i) ))
      } else {
        mat[a+1,a+1]=mat[a+1,a+1]-0.25
        b=bitwXor(a,2^i+2^j)
        mat[a+1,b+1]=mat[a+1,b+1]+0.5
#        print(paste( "mat(a,a)=", mat[a+1,a+1], " b=", b, " mat(a,b)=", mat[a+1,b+1] , btest(a,i), btest(a,j)  ))
      }
    }
  }

  return(mat)
}

mat = hamiltonian(nn, nst, mat)


