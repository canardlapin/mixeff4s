d <- data.frame(f=factor(c('a','a','b','b')), h=factor(c('u','v','u','v')), x=1:4)
witnesses <- list(
  list(~0+f+h, matrix(c(1,1,0,0, 0,0,1,1, 0,1,0,1),4)),
  list(~0+f:h, matrix(c(1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1),4)),
  list(~h+f:h, matrix(c(1,1,1,1, 0,1,0,1, 0,0,1,0, 0,0,0,1),4)),
  list(~0+x:f, matrix(c(1,2,0,0, 0,0,3,4),4))
)
print(R.version.string)
for (w in witnesses) {
  actual <- model.matrix(w[[1]], d)
  stopifnot(isTRUE(all.equal(unname(actual), w[[2]], check.attributes=FALSE)))
  print(w[[1]])
  print(actual)
}
