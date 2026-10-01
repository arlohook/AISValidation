library(fda)
library(mvtnorm)

#take base data and generate from the covariance matrix
t = seq(0,1, 0.01)
X = eval.fd(t, fd(gait[,,2], create.bspline.basis(c(0,1), 20, 4)))
matplot(X, type = 'l')

n = ncol(X)
Xmu = apply(X, 1, mean)
Xc = apply(X,2, function(x){x-Xmu})


Ci = chol(tcrossprod(Xc)/(n-1)+ diag(1e-10, 101))

set.seed(2026)

xi = rmvnorm(100, mean = rep(0,101))
Xchat = t(xi%*%Ci)
Xhat = apply(Xchat, 2, function(x){x+Xmu})



matplot(Xhat, type = 'l')

# make trial reference
ref_df = data.frame(Participant = rep(1:10, each = 10),
                    Trial = 1:100)

# make participant specific deviations
sigma = 1
ell = 0.3

# pairwise squared distances
D2 <- outer(t, t, function(a, b) (a - b)^2)
# Gaussian (squared-exponential) kernel
C <- sigma^2 * exp(-D2 / (2 * ell^2))

fields::image.plot(C)

jdif = t(rmvnorm(10, sigma = C))

matplot(jdif, type = 'l')

# make global deviation
alpha = 3*cos(2.2*t*pi+15)
plot(alpha)+
  abline(h = 0)

jdmat = do.call(cbind, lapply(rep(1:10, each = 10), function(j){jdif[,j]}))


# make smooth noise
sigma = 1
ell = 0.2
Cep <- sigma^2 * exp(-D2 / (2 * ell^2))

fields::image.plot(Cep)

noise = t(rmvnorm(100, sigma = Cep))

# combine to make second system
Xhat2 = apply(Xhat, 2, function(i){i+alpha}) +jdmat + noise

par(mfrow = c(1,3))
matplot(Xhat, main = "System 1", type = 'l')
matplot(Xhat2, main = "System 2", type = 'l')
matplot(Xhat-Xhat2, typ = 'l', main = "Difference")


# combine and save
out = list("Reference" = ref_df,
           "System1" = Xhat,
           "System2" = Xhat2)

saveRDS(out, "Demo Data.rds")





