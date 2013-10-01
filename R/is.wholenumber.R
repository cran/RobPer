is.wholenumber <-
function(x, tol = .Machine$double.eps^0.5)  abs(abs(x) - round(x)) < tol
