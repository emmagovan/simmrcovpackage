# #This is just to test stuff in the package
# 
# a = c(1:5)
# b = c(6:10)
# c = c(11:15)
# d = c(16:20)
# 
# 
# #This doesn't work because it just creates a vector out of a and b
# formula = c(a,b) ~ c * d
# 
# model.matrix(formula)
# 
# #We need it to be a matrix if they are our two isotopes
# #e.g.
# 
# abmat = matrix(c(a,b), ncol = 2)
# 
# formula = abmat ~ c * d
# 
# model.matrix(formula)
