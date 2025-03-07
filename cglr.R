# R code for the cglr package by Marco Geraci, Sapienza University of Rome

#########################################################
### Required packages
#########################################################

pckgs <- c("nlme", "Rcpp", "RcppEigen", "RcppArmadillo", "SphericalCubature", "nlmm", "cubature", "mvtnorm", "numDeriv", "corpcor", "boot")

# Install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Apply function to each package
sapply(pckgs, install_if_missing)

# Load all packages
lapply(pckgs, library, character.only = TRUE)

#########################################################
### Spare functions
#########################################################

cor2cov <- function(R, sd) {
  R * tcrossprod(sd)
}

clock2rad <- function(h, m, s = 0){
(h + m/60 + s/3600)*pi/12
}

rad2clock <- function(x, decimal = FALSE, counter = FALSE){
x <- ifelse(counter, 2*pi + x, x)
val <- x*12/pi
h <- floor(val)
m <- floor((val-h)*60)
s <- floor(((val-h)*60-m)*60)
out <- list(h = h, m = m, s = s)
if(decimal) out <- h+m/60+s/3600
return(out)
}

#########################################################
### Cylindrical normal regression
#########################################################

cpnrControl <- function(optimizer = "nlminb", method = "Nelder-Mead", maxit = 500, verbose = FALSE) 
{
    list(optimizer = optimizer, method = method, maxit = maxit, verbose = verbose)
}

loglik_cpnr <- function(eta, w, theta, r, mmx, conditional){

# dimensions
n <- length(theta)
K <- ncol(w)
P1 <- ncol(mmx)*(K+2)
P2 <- (K+1) + (K+1)*(K+2)/2 # phi and Phi, equation (8)
if(conditional) P2 <- P2 + 1 # if conditional, variance is not constrained

# get parameters
b <- eta[1:P1]
xb <- mmx %*% matrix(b, ncol = K + 2)
gamma <- eta[(P1+1):(P1+P2)]
sigma <- buildSigma(gamma, K = K, conditional = conditional)

if(conditional){
	y <- cbind(w, r*cos(theta), r*sin(theta))
	val <- C_loglik_cond_cpnr(y = as.matrix(y), m = as.matrix(xb), sigma = as.matrix(sigma))
} else {
	val <- C_loglik_proj_cpnr(w = as.matrix(w), theta = as.vector(theta), m = as.matrix(xb), sigma = as.matrix(sigma))
}

# negative log-likelihood
return(-val)
}

summary.cpnr <- function(object, alpha = 0.05, ...){

ynn <- if(!is.null(object$response)) c(object$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")

	g <- function(x){
		exp(x)
	}

K <- object$K
b <- object$b
B <- matrix(b, ncol = K + 2)
P <- nrow(B)
P1 <- object$dims[["P1"]]
P2 <- object$dims[["P2"]]

V <- vcov(object)

SE_b <- sqrt(diag(V)[1:P1])
lower_b <- b - SE_b * qnorm(1 - alpha/2, 0, 1)
upper_b <- b + SE_b * qnorm(1 - alpha/2, 0, 1)

tTable_b <- list()
for(j in 1:(K+2)){
	tTable_b[[j]] <- data.frame(B[,j], matrix(SE_b, nrow = P)[,j], matrix(lower_b, nrow = P)[,j], matrix(upper_b, nrow = P)[,j])
	names(tTable_b[[j]]) <- c("Estimate", "Std.Err", "Lower", "Upper")
	rownames(tTable_b[[j]]) <- object$terms
}
names(tTable_b) <- ynn
object$tTable_b <- tTable_b
class(object) <- "summary.cpnr"
return(object)
}

print.summary.cpnr <- function(x, digits = max(3, getOption("digits") - 3), ...){

K <- x$K
Sigma <- buildSigma(x$gamma, K = K, conditional = x$conditional)

ynn <- if(!is.null(x$response)) c(x$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")

colnames(Sigma) <- rownames(Sigma) <- ynn

cat("Cylindrical normal model ")
if(x$conditional) cat("(conditional on circle radius)", "\n")
if(!x$conditional) cat("(projected on circle)", "\n")
cat("Log-likelihood:\n")
print.default(format(x$loglik, digits = digits), print.gap = 2,	quote = FALSE)
cat("Fixed effects:\n")
print(x$tTable_b, ...)
cat("\n")
cat("Scale matrix of the error:\n")
print.default(format(as.matrix(Sigma), digits = digits), quote = FALSE)
cat("\n")
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat("\n")
invisible(x)
}

cpnr <- function(formula, data, r = NULL, conditional = FALSE, control = NULL, start = NULL, fit = TRUE){

Call <- match.call()

# parse formula
w <- model.matrix(as.formula(paste("~ 0 +", deparse(formula[[2]][[2]]))), data)
theta <- model.matrix(as.formula(paste("~ 0 +", deparse(formula[[2]][[3]]))), data)
mmx <- model.matrix(as.formula(paste("~ ", deparse(formula[[3]]))), data)

# dimensions
n <- length(theta)
K <- ncol(w)
P1 <- ncol(mmx)*(K+2)
P2 <- (K+1) + (K+1)*(K+2)/2 # phi and Phi, equation (8)
if(conditional) P2 <- P2 + 1 # if conditional, variance is not constrained

# optimization control parameters
if (is.null(names(control))) 
	control <- cpnrControl()
else {
	control_default <- cpnrControl()
	control_names <- intersect(names(control), names(control_default))
	control_default[control_names] <- control[control_names]
	control <- control_default
}

# initialize
if(is.null(start)){
	y <- cbind(w, cos(theta), sin(theta))
	dimnames(y) <- NULL
	B <- lm(y ~ mmx - 1)$coef
	b <- as.vector(B)
	tmp <- var(y)
	tmp <- tmp/tmp[K+2,K+2]
	gamma <- extractSigma(tmp, K = K, conditional = conditional)
	start <- eta_0 <- c(b, gamma)
} else {
	dim_eta <- P1+P2
	if(length(start) != dim_eta) stop("start is of length ", length(start), ". Should be ", dim_eta)
	eta_0 <- start
}

if(is.null(r)) r <- rep(1, n)

FIT_ARGS <- list(eta = eta_0, w = w, theta = theta, r = r, mmx = mmx, conditional = conditional)

if(!fit) return(FIT_ARGS)

if(control$optimizer == "nlminb"){
	fit <- do.call(nlminb, args = c(list(objective = loglik_cpnr, start = FIT_ARGS$eta, control = list(trace = as.numeric(control$verbose))), FIT_ARGS[-c(match(c("eta"), names(FIT_ARGS)))]))
	cat(fit$message, "\n")
	fit$loglik <- -fit$objective
}

if(control$optimizer == "optim"){
	fit <- do.call(optim, args = c(list(fn = loglik_cpnr, par = FIT_ARGS$eta, method = control$method, control = list(maxit = control$maxit, trace = as.numeric(control$verbose))), FIT_ARGS[-c(match(c("eta"), names(FIT_ARGS)))]))
	if(fit$convergence == 0) cat("Convergence successful", "\n")
	if(fit$convergence == 1) cat("Reached maximum number of iterations without convergence", "\n")
	fit$loglik <- -fit$value
}

names(fit$par) <- NULL
fit$call <- Call
fit$K <- K
fit$b <- fit$par[1:P1]
fit$gamma <- fit$par[(P1+1):(P1+P2)]
fit$conditional <- conditional
fit$control <- control
fit$w <- w
fit$theta <- theta
fit$start <- start
fit$r <- r
fit$mmx <- mmx
fit$terms <- colnames(mmx)
fit$response <- colnames(w)
fit$nobs <- n
fit$dims <- c(P1 = P1, P2 = P2)
fit$formula <- formula

class(fit) <- "cpnr"
return(fit)
}

vcov.cpnr <- function(object, ...){

H <- numDeriv::hessian(func = loglik_cpnr, x = object$par, w = object$w, theta = object$theta, r = object$r, mmx = object$mmx, conditional = object$conditional)
H[is.nan(H)] <- 1e-6

val <- solve(H)
if(!corpcor::is.positive.definite(val)) val <- corpcor::make.positive.definite(val)
val
}

predict.cpnr <- function(object, obs = -1){

	K <- object$K
	ynn <- if(!is.null(object$response)) c(object$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")
	
	n <- length(object$theta)
	B <- coef(object)
	
	if(obs < 0){
		val <- object$mmx %*% B
		colnames(val) <- ynn
	}
	if(obs == 0){
		xm <- apply(object$mmx, 2, mean)
		val <- as.numeric(matrix(xm, nrow = 1) %*% B)
		names(val) <- ynn
	}
	if(obs > 0){
		if(length(obs) != 1) stop("One observation at a time or set 'obs = -1' for all observations")
		if(obs > n) stop("Observation must be ", n, " at most")
		val <- as.numeric(object$mmx[obs,,drop=FALSE] %*% B)
		names(val) <- ynn
	}
	
	return(val)
}

coef.cpnr <- function(object){

	K <- object$K
	ynn <- if(!is.null(object$response)) c(object$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")
	
	val <- matrix(object$b, ncol = K + 2)
	colnames(val) <- ynn
	rownames(val) <- object$terms

	return(val)
}

print.cpnr <- function(x, digits = max(3, getOption("digits") - 3), ...){

K <- x$K
B <- matrix(x$b, ncol = K + 2)
Sigma <- buildSigma(x$gamma, K = K, conditional = x$conditional)

ynn <- if(!is.null(x$response)) c(x$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")
xnn <- x$terms

colnames(B) <- ynn
rownames(B) <- xnn
colnames(Sigma) <- rownames(Sigma) <- ynn

cat("Cylindrical normal model ")
if(x$conditional) cat("(conditional on circle radius)", "\n")
if(!x$conditional) cat("(projected on circle)", "\n")
cat("Log-likelihood:\n")
print.default(format(x$loglik, digits = digits), print.gap = 2,	quote = FALSE)
cat("Fixed effects:\n")
print.default(format(B, digits = digits), print.gap = 2, quote = FALSE)
cat("\n")
cat("Scale matrix of the error:\n")
print.default(format(as.matrix(Sigma), digits = digits), quote = FALSE)
cat("\n")
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat("\n")
invisible(x)
}

getVarCov.cpnr <- function(obj){

K <- obj$K
n <- length(obj$theta)
ynn <- if(!is.null(obj$response)) c(obj$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")

sigma <- buildSigma(obj$gamma, K = K, conditional = obj$conditional)
colnames(sigma) <- rownames(sigma) <- ynn

return(sigma)
}

corlc.cpnr <- function(object){

K <- object$K
n <- length(object$theta)
ynn <- if(!is.null(object$response)) object$response else paste("Linear outcome", 1:K)

sigma <- getVarCov(object)
ans <- rep(NA, K)
Rho <- cov2cor(sigma)
for(j in 1:K){
	rxc <- Rho[j,(K+1)]
	rxs <- Rho[j,(K+2)]
	rcs <- Rho[(K+1),(K+2)]
	ans[j] <- ((rxc*rxc)+(rxs*rxs)-(2*rxc*rxs*rcs))/(1-rcs*rcs)
}
names(ans) <- ynn

return(ans)
}

residuals.cpnr <- function(object, ...){

K <- object$K

ynn <- if(!is.null(object$response)) c(object$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")
	
yhat <- predict(object, obs = -1)
y <- cbind(object$w, cos(object$theta), sin(object$theta))
ans <- y - yhat
colnames(ans) <- ynn

return(ans)
}

logLik.cpnr <- function(object, ...){
val <- object$loglik
attr(val, "df") <- sum(object$dims)
attr(val, "rdf") <- object$nobs - sum(object$dims)
return(val)
}

AIC.cpnr <- function(object, ..., k = 2){
p <- sum(object$dims)
val <- 2*p - 2*object$loglik
return(val)
}

# Get marginal densities from cpnr object

dtheta.cpnr <- function(object, data = NULL, r = NULL){

K <- object$K

if(is.null(data)){
	theta <- object$theta
	mmx <- object$mmx
} else {
	theta <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[3]]))), data)
	mmx <- model.matrix(as.formula(paste("~ ", deparse(object$formula[[3]]))), data)
}

n <- length(theta)
if(object$conditional & is.null(r)){
	r <- object$r
}

if(object$conditional){
	if(length(r) != n) warning("Length of r (", length(r), ") is not the same as other variables (", n, ")")
}

m <- mmx %*% coef(object)
m <- m[,-c(1:K),drop=FALSE]
sigma <- getVarCov(object)
sigmas <- sigma[-c(1:K),-c(1:K)]

if(object$conditional){
	#tmp <- try(zzz(as.vector(theta), as.matrix(m), as.matrix(sigmas), as.vector(r)), silent = TRUE)
} else {
	val <- try(C_dpn_vv(theta = as.vector(theta), m = as.matrix(m), sigma = as.matrix(sigmas), logd = as.logical(FALSE)), silent = TRUE)
}
return(val)
}

dw.cpnr <- function(object, variable, data = NULL){

K <- object$K

if(is.null(data)){
	w <- object$w
	mmx <- object$mmx
} else {
	w <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[2]]))), data)
	mmx <- model.matrix(as.formula(paste("~ ", deparse(object$formula[[3]]))), data)
}
n <- nrow(w)
if(ncol(w) != K) stop("w must have K =", K, "columns")

m <- mmx %*% coef(object)
m <- m[,c(1:K),drop=FALSE]
sigma <- getVarCov(object)
sigmaw <- sigma[1:K,1:K,drop=FALSE]

if(variable == "theta") stop("Use dtheta.cpnr or dens.cpnr(..., variable = 'theta')")
sel <- match(variable, object$response)
if(is.na(sel)) stop("variable invalid")

val <- rep(NA, n)
for(i in 1:n){
	z <- as.vector(w[i,] - m[i,])
	val[i] <- try(mvtnorm::dmvnorm(z[sel], sigma = sigmaw[sel,sel,drop=FALSE], log = FALSE), silent = TRUE)
}
return(val)
}

dens.cpnr <- function(object, variable = "theta", data = NULL){

K <- object$K

g <- function(xnsel, xsel, sel, nsel, K, m, sigma, muw, alpha){
	x <- rep(0, K+1)
	x[sel] <- xsel
	x[nsel] <- xnsel
	w <- x[1:K]
	theta <- x[K+1]
	ans <- C_loglik_proj_cpnr(w = matrix(w, nrow = 1), theta = as.vector(theta), m = matrix(m, nrow = 1), sigma = as.matrix(sigma))
	exp(ans)
}

if(is.null(data)){
	w <- object$w
	theta <- object$theta
	mmx <- object$mmx
} else {
	w <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[2]]))), data)
	theta <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[3]]))), data)
	mmx <- model.matrix(as.formula(paste("~ ", deparse(object$formula[[3]]))), data)
}
n <- nrow(w)
if(ncol(w) != K) stop("w must have K =", K, "columns")

m <- mmx %*% coef(object)
sigma <- getVarCov(object)

sel <- match(variable, c(object$response, "theta"))
if(is.na(sel)) stop("variable invalid")
nsel <- setdiff(1:(K+1), sel)

y <- cbind(w, theta)
low <- apply(y, 2, min)
up <- apply(y, 2, max)
low[K+1] <- 0
up[K+1] <- 2*pi
low <- low[nsel]
up <- up[nsel]

val <- rep(NA, n)
for(i in 1:n){
	val[i] <- try(cubature::pcubature(f = g, lowerLimit = low, upperLimit = up, xsel = y[i,sel], sel = sel, nsel = nsel, K = K, m = m[i,], sigma = sigma, tol = 1e-03, vectorInterface = FALSE)$integral, silent = TRUE)
}
return(val)
}

#########################################################
### Cylindrical GL regression
#########################################################

cglrControl <- function(optimizer = "nlminb", method = "Nelder-Mead", maxit = 500, Hv = 7, multistart = TRUE, grid = c(0.1, 1, 2), alpha = c(1, 1), alpha.index = 9, verbose = FALSE) 
{
    if (length(alpha) != 2) 
        stop("Provide starting values for alpha")
    if (!alpha.index %in% c(0, 1, 2, 9)) 
        stop("alpha.index is one of c(0,1,2,9)")
    if (any(alpha < 0)) 
        stop("values for alpha must be greater than 0")
    if (any(grid < 0)) 
        stop("values for alpha.grid must be greater than 0")
    list(optimizer = optimizer, method = method, maxit = maxit, Hv = Hv, multistart = multistart, grid = grid, alpha = alpha, alpha.index = alpha.index, verbose = verbose)
}

buildSigma <- function(x, K = NULL, conditional = FALSE){

if(conditional){
	val <- as.matrix(nlme::pdSymm(value = x))
} else {
	phi <- matrix(x[1:(K+1)])
	Phi <- as.matrix(nlme::pdSymm(value = x[-c(1:(K+1))]))
	up <- cbind(Phi + phi%*%t(phi), phi)
	low <- c(as.numeric(phi), 1)
	val <- rbind(up, low)
	dimnames(val) <- NULL
}

if(!corpcor::is.positive.definite(val)){
warning("Matrix is not pd")
}

return(val)
}

extractSigma <- function(x, K = NULL, conditional = FALSE){

if(conditional){
	val <- coef(nlme::pdSymm(x))
} else {
	Phi <- coef(nlme::pdSymm(x[1:(K+1),1:(K+1)]))
	phi <- x[K+2,1:(K+1)]
	val <- c(phi, Phi)
}

return(val)
}

loglik_cglr <- function(eta, w, theta, r, mmx, mmu, mma, Hv, conditional){

# dimensions
n <- length(theta)
K <- ncol(w)
P1 <- ncol(mmx)*(K+2)
P2 <- (K+1) + (K+1)*(K+2)/2 # phi and Phi, equation (8)
if(conditional) P2 <- P2 + 1 # if conditional, variance is not constrained
P3 <- ncol(mmu)*K
P4 <- ncol(mma)

# get parameters
b <- eta[1:P1]
xb <- mmx %*% matrix(b, ncol = K + 2)
gamma <- eta[(P1+1):(P1+P2)]
Sigma <- buildSigma(gamma, K = K, conditional = conditional)
muw <- eta[(P1+P2+1):(P1+P2+P3)]
xmu <- mmu %*% matrix(muw, ncol = K)
tau <- eta[(P1+P2+P3+1):(P1+P2+P3+P4)]
if(all(colnames(mma) == "(Intercept)")){
	alpha <- exp(tau[1])
} else {
	xtau <- mma %*% tau
	alpha <- exp(xtau)
	alpha <- pmin(alpha, 1e10)
}

if(conditional){
	y <- cbind(w, r*cos(theta), r*sin(theta))
	val <- C_loglik_cond_cglr(y = as.matrix(y), xb = as.matrix(xb), sigma = as.matrix(Sigma), mu = as.matrix(cbind(xmu, 0, 0)), alpha = as.vector(alpha), Hv = as.integer(Hv))
} else {
	val <- C_loglik_proj_cglr(w = as.matrix(w), theta = as.vector(theta), xb = as.matrix(xb), sigma = as.matrix(Sigma), muw = as.matrix(xmu), alpha = as.vector(alpha), Hv = as.integer(Hv))
}

# negative integrated log-likelihood
return(-sum(val))
}

summary.cglr <- function(object, alpha = 0.05, ...){

ynn <- if(!is.null(object$response)) c(object$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")

	g <- function(x){
		exp(x)
	}

K <- object$K
b <- object$b
B <- matrix(b, ncol = K + 2)
P <- nrow(B)
P1 <- object$dims[["P1"]]
P2 <- object$dims[["P2"]]
P3 <- object$dims[["P3"]]
P4 <- object$dims[["P4"]]

tau <- object$tau
alphap <- exp(tau)
Sigma <- buildSigma(object$gamma, K = K, conditional = object$conditional)
muw <- as.vector(object$muw)
mu <- cbind(matrix(object$muw, ncol = K), 0, 0)

V <- vcov(object)

SE_b <- sqrt(diag(V)[1:P1])
lower_b <- b - SE_b * qnorm(1 - alpha/2, 0, 1)
upper_b <- b + SE_b * qnorm(1 - alpha/2, 0, 1)

SE_tau <- sqrt(diag(V)[(P1+P2+P3+1):(P1+P2+P3+P4)])
SE_alpha <- sqrt(SE_tau^2*g(tau)^2)
lower_tau <- tau - SE_tau * qnorm(1 - alpha/2, 0, 1)
upper_tau <- tau + SE_tau * qnorm(1 - alpha/2, 0, 1)
lower_alpha <- exp(lower_tau)
upper_alpha <- exp(upper_tau)

SE_mu <- sqrt(diag(V)[(P1+P2+1):(P1+P2+P3)])
lower_mu <- muw - SE_mu * qnorm(1 - alpha/2, 0, 1)
upper_mu <- muw + SE_mu * qnorm(1 - alpha/2, 0, 1)

m1 <- cbind(matrix(SE_mu, ncol = K), NA, NA)
m2 <- cbind(matrix(lower_mu, ncol = K), NA, NA)
m3 <- cbind(matrix(upper_mu, ncol = K), NA, NA)

tTable_alpha <- data.frame(alphap, c(SE_alpha), c(lower_alpha), c(upper_alpha))
names(tTable_alpha) <- c("Estimate", "Std.Err", "Lower", "Upper")
rownames(tTable_alpha) <- object$terms[["mma"]]

tTable_b <- list()
tTable_mu <- list()
for(j in 1:(K+2)){
	tTable_b[[j]] <- data.frame(B[,j], matrix(SE_b, nrow = P)[,j], matrix(lower_b, nrow = P)[,j], matrix(upper_b, nrow = P)[,j])
	names(tTable_b[[j]]) <- c("Estimate", "Std.Err", "Lower", "Upper")
	rownames(tTable_b[[j]]) <- object$terms[["mmx"]]
	
	tTable_mu[[j]] <- data.frame(mu[,j], matrix(m1, nrow = ncol(object$mmu))[,j], matrix(m2, nrow = ncol(object$mmu))[,j], matrix(m3, nrow = ncol(object$mmu))[,j])
	names(tTable_mu[[j]]) <- c("Estimate", "Std.Err", "Lower", "Upper")
	rownames(tTable_mu[[j]]) <- object$terms[["mmu"]]
}
names(tTable_b) <- ynn
names(tTable_mu) <- ynn
object$tTable_b <- tTable_b
object$tTable_mu <- tTable_mu
object$tTable_alpha <- tTable_alpha
class(object) <- "summary.cglr"
return(object)
}

print.summary.cglr <- function(x, digits = max(3, getOption("digits") - 3), ...){

K <- x$K
Sigma <- buildSigma(x$gamma, K = K, conditional = x$conditional)

ynn <- if(!is.null(x$response)) c(x$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")

colnames(Sigma) <- rownames(Sigma) <- ynn

cat("Cylindrical Generalized Laplace model ")
if(x$conditional) cat("(conditional on circle radius)", "\n")
if(!x$conditional) cat("(projected on circle)", "\n")
cat("Log-likelihood:\n")
print.default(format(x$loglik, digits = digits), print.gap = 2,	quote = FALSE)
cat("Alpha:\n")
print(x$tTable_alpha, ...)
cat("\n")
cat("Fixed effects:\n")
print(x$tTable_b, ...)
cat("\n")
cat("Asymmetry:\n")
print(x$tTable_mu, ...)
cat("\n")
cat("Scale matrix of the error:\n")
print.default(format(as.matrix(Sigma), digits = digits), quote = FALSE)
cat("\n")
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat("\n")
invisible(x)
}

vcov.cglr <- function(object, ...){

H <- numDeriv::hessian(func = loglik_cglr, x = object$par, w = object$w, theta = object$theta, r = object$r, mmx = object$mmx, mmu = object$mmu, mma = object$mma, Hv = object$control$Hv, conditional = object$conditional)
H[is.nan(H)] <- 1e-6

val <- solve(H)
if(!corpcor::is.positive.definite(val)) val <- corpcor::make.positive.definite(val)
val
}

logLik.cglr <- function(object, ...){
val <- object$loglik
attr(val, "df") <- sum(object$dims)
attr(val, "rdf") <- object$nobs - sum(object$dims)
return(val)
}

AIC.cglr <- function(object, ..., k = 2){
p <- sum(object$dims)
val <- 2*p - 2*object$loglik
return(val)
}

cglr <- function(formula, mu.formula = ~ 1, alpha.formula = ~ 1, data, r = NULL, conditional = FALSE, control = NULL, start = NULL, fit = TRUE){

Call <- match.call()

# parse formula
w <- model.matrix(as.formula(paste("~ 0 +", deparse(formula[[2]][[2]]))), data)
theta <- model.matrix(as.formula(paste("~ 0 +", deparse(formula[[2]][[3]]))), data)
mmx <- model.matrix(as.formula(paste("~ ", deparse(formula[[3]]))), data)

# mu formula
mmu <- model.matrix(mu.formula, data)

# alpha formula
mma <- model.matrix(alpha.formula, data)

# dimensions
n <- length(theta)
K <- ncol(w)
P1 <- ncol(mmx)*(K+2)
P2 <- (K+1) + (K+1)*(K+2)/2 # phi and Phi, equation (8)
if(conditional) P2 <- P2 + 1 # if conditional, variance is not constrained
P3 <- ncol(mmu)*K
P4 <- ncol(mma)

# optimization control parameters
if (is.null(names(control))) 
	control <- cglrControl()
else {
	control_default <- cglrControl()
	control_names <- intersect(names(control), names(control_default))
	control_default[control_names] <- control[control_names]
	control <- control_default
}

# initialize
if(is.null(start)){
	y <- cbind(w, cos(theta), sin(theta))
	dimnames(y) <- NULL
	B <- lm(y ~ mmx - 1)$coef
	b <- as.vector(B)
	tmp <- var(y)
	tmp <- tmp/tmp[K+2,K+2]
	gamma <- extractSigma(tmp, K = K, conditional = conditional)
	muw <- rep(0, P3)
	mu <- c(muw, 0, 0)
	tau <- rep(0, P4)
	start <- eta_0 <- c(b, gamma, muw, tau)
} else {
	dim_eta <- P1+P2+P3+P4
	if(length(start) != dim_eta) stop("start is of length ", length(start), ". Should be ", dim_eta)
	eta_0 <- start
}

if(is.null(r)) r <- rep(1, n)


FIT_ARGS <- list(eta = eta_0, w = w, theta = theta, r = r, mmx = mmx, mmu = mmu, mma = mma, Hv = control$Hv, conditional = conditional)

if(!fit) return(FIT_ARGS)

if(control$optimizer == "nlminb"){
	fit <- do.call(nlminb, args = c(list(objective = loglik_cglr, start = FIT_ARGS$eta, control = list(trace = as.numeric(control$verbose))), FIT_ARGS[-c(match(c("eta"), names(FIT_ARGS)))]))
	cat(fit$message, "\n")
	fit$loglik <- -fit$objective
}

if(control$optimizer == "optim"){
	fit <- do.call(optim, args = c(list(fn = loglik_cglr, par = FIT_ARGS$eta, method = control$method, control = list(maxit = control$maxit, trace = as.numeric(control$verbose))), FIT_ARGS[-c(match(c("eta"), names(FIT_ARGS)))]))
	if(fit$convergence == 0) cat("Convergence successful", "\n")
	if(fit$convergence == 1) cat("Reached maximum number of iterations without convergence", "\n")
	fit$loglik <- -fit$value
}

names(fit$par) <- NULL
fit$call <- Call
fit$K <- K
fit$b <- fit$par[1:P1]
fit$gamma <- fit$par[(P1+1):(P1+P2)]
fit$muw <- fit$par[(P1+P2+1):(P1+P2+P3)]
fit$tau <- fit$par[(P1+P2+P3+1):(P1+P2+P3+P4)]
fit$conditional <- conditional
fit$control <- control
fit$w <- w
fit$theta <- theta
fit$start <- start
fit$r <- r
fit$mmx <- mmx
fit$mmu <- mmu
fit$mma <- mma
fit$terms <- list(mmx = colnames(mmx), mmu = colnames(mmu), mma = colnames(mma))
fit$response <- colnames(w)
fit$nobs <- n
fit$dims <- c(P1 = P1, P2 = P2, P3 = P3, P4 = P4)
fit$formula <- formula
fit$mu.formula <- mu.formula
fit$alpha.formula <- alpha.formula

class(fit) <- "cglr"
return(fit)
} #

getVarCov.cglr <- function(obj, obs = 0, scale = FALSE){

K <- obj$K
n <- length(obj$theta)
ynn <- if(!is.null(obj$response)) c(obj$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")

mu <- cbind(matrix(obj$muw, ncol = K), 0, 0)
Sigma <- buildSigma(obj$gamma, K = K, conditional = obj$conditional)
colnames(Sigma) <- rownames(Sigma) <- ynn

if(obs < 0){
	val <- array(NA, dim = c(K+2, K+2, n), dimnames = list(ynn, ynn, 1:n))
	xmu <- obj$mmu %*% mu
	alpha <- exp(obj$mma %*% obj$tau)
	for(i in 1:n){
		mm <- crossprod(xmu[i,,drop=FALSE])
		val[,,i] <- if(scale) Sigma else (Sigma + mm)*alpha[i]
	}
}
if(obs == 0){
	xm <- apply(obj$mma, 2, mean)
	alpha <- exp(sum(xm * obj$tau))
	xm <- apply(obj$mmu, 2, mean)
	xmu <- as.numeric(matrix(xm, nrow = 1) %*% mu)
	mm <- tcrossprod(xmu)
	val <- if(scale) Sigma else (Sigma + mm)*alpha
	colnames(val) <- rownames(val) <- ynn
}
if(obs > 0){
	if(any(obs > n)) stop("Observation must be ", n, " at most")
	alpha <- exp(sum(obj$mma[obs,,drop=TRUE] * obj$tau))
	xmu <- as.numeric(obj$mmu[obs,,drop=FALSE] %*% mu)
	mm <- tcrossprod(xmu)
	val <- if(scale) Sigma else (Sigma + mm)*alpha
	colnames(val) <- rownames(val) <- ynn
}

return(val)
}

corlc.cglr <- function(object, obs = 0){

K <- object$K
n <- length(object$theta)
ynn <- if(!is.null(object$response)) object$response else paste("Linear outcome", 1:K)

Sigma <- getVarCov(object, obs = obs, scale = FALSE)
if(obs < 0){
	ans <- array(NA, dim = c(n, K), dimnames = list(1:n, ynn))
	for(i in 1:n){
	Rho <- cov2cor(Sigma[,,i])
		for(j in 1:K){
			rxc <- Rho[j,(K+1)]
			rxs <- Rho[j,(K+2)]
			rcs <- Rho[(K+1),(K+2)]
			ans[i,j] <- ((rxc*rxc)+(rxs*rxs)-(2*rxc*rxs*rcs))/(1-rcs*rcs)
		}
	}
} else {
	ans <- rep(NA, K)
	Rho <- cov2cor(Sigma)
	for(j in 1:K){
		rxc <- Rho[j,(K+1)]
		rxs <- Rho[j,(K+2)]
		rcs <- Rho[(K+1),(K+2)]
		ans[j] <- ((rxc*rxc)+(rxs*rxs)-(2*rxc*rxs*rcs))/(1-rcs*rcs)
	}
	names(ans) <- ynn
}

return(ans)
}

residuals.cglr <- function(object, ...){

K <- object$K

ynn <- if(!is.null(object$response)) c(object$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")
	
yhat <- predict(object, obs = -1, location = TRUE)
y <- cbind(object$w, cos(object$theta), sin(object$theta))
ans <- y - yhat
colnames(ans) <- ynn

return(ans)
}

predict.cglr <- function(object, obs = -1, location = FALSE){

	K <- object$K
	ynn <- if(!is.null(object$response)) c(object$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")
	
	n <- length(object$theta)
	B <- coef(object, "beta")
	mu <- coef(object, "mu")
	
	if(obs < 0){
		delta <- object$mmx %*% B
		xmu <- object$mmu %*% mu
		alpha <- exp(object$mma %*% object$tau)
		val <- if(location) delta else delta + sweep(xmu, 1, alpha, FUN = "*")
		colnames(val) <- ynn
	}
	if(obs == 0){
		xm <- apply(object$mmx, 2, mean)
		delta <- as.numeric(matrix(xm, nrow = 1) %*% B)
		xm <- apply(object$mmu, 2, mean)
		xmu <- as.numeric(matrix(xm, nrow = 1) %*% mu)
		xm <- apply(object$mma, 2, mean)
		alpha <- exp(sum(xm * object$tau))
		val <- if(location) delta else delta + alpha*xmu
		names(val) <- ynn
	}
	if(obs > 0){
		if(length(obs) != 1) stop("One observation at a time or set 'obs = -1' for all observations")
		if(obs > n) stop("Observation must be ", n, " at most")
		delta <- as.numeric(object$mmx[obs,,drop=FALSE] %*% B)
		xmu <- as.numeric(object$mmu[obs,,drop=FALSE] %*% mu)
		alpha <- exp(sum(object$mma[obs,,drop=TRUE] * object$tau))
		val <- if(location) delta else delta + alpha*xmu
		names(val) <- ynn
	}
	
	return(val)
}

coef.cglr <- function(object, which = "beta"){

	K <- object$K
	ynn <- if(!is.null(object$response)) c(object$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")
	
	if(which == "beta"){
		val <- matrix(object$b, ncol = K + 2)
		colnames(val) <- ynn
		rownames(val) <- object$terms[["mmx"]]
	}
	if(which == "mu"){
		val <- cbind(matrix(object$muw, ncol = K), 0, 0)
		colnames(val) <- ynn
		rownames(val) <- object$terms[["mmu"]]
	}
	if(which == "tau"){
		val <- object$tau
		names(val) <- object$terms[["mma"]]
	}
	
	
	return(val)
}

boot.cglr <- function(object, data, R = 100, seed = 1){

	f <- function(data, indices, object) {
		K <- object$K
		fitted_values <- predict(object, obs = -1, location = TRUE)
		res <- residuals(object)
		y <- fitted_values + res[indices]
		
		FIT_ARGS <- object[c("r", "mmx", "mmu", "mma", "conditional")]
		FIT_ARGS$w <- y[,c(1:K),drop=FALSE]
		FIT_ARGS$theta <- atan2(y[,K+2],y[,K+1])
		FIT_ARGS$Hv <- object$control$Hv
		fit <- do.call(optim, args = c(list(fn = loglik_cglr, par = object$par, method = object$control$method, control = list(maxit = object$control$maxit)), FIT_ARGS))
		return(fit$par)
	}

set.seed(seed)
ans <- boot::boot(data = data, statistic = f, R = R, object = object)

return(ans)

}

print.cglr <- function(x, digits = max(3, getOption("digits") - 3), ...){

K <- x$K
B <- matrix(x$b, ncol = K + 2)
tau <- x$tau
alpha <- exp(tau)
Sigma <- buildSigma(x$gamma, K = K, conditional = x$conditional)
mu <- cbind(matrix(x$muw, ncol = K), 0, 0)

ynn <- if(!is.null(x$response)) c(x$response, "cos(theta)", "sin(theta)") else c(paste("Linear outcome", 1:K), "cos(theta)", "sin(theta)")
xnn <- x$terms[['mmx']]

colnames(B) <- ynn
rownames(B) <- x$terms[['mmx']]
colnames(Sigma) <- rownames(Sigma) <- ynn
colnames(mu) <- ynn
rownames(mu) <- x$terms[['mmu']]
names(tau) <- names(alpha) <- x$terms[['mma']]

cat("Cylindrical Generalized Laplace model ")
if(x$conditional) cat("(conditional on circle radius)", "\n")
if(!x$conditional) cat("(projected on circle)", "\n")
cat("Log-likelihood:\n")
print.default(format(x$loglik, digits = digits), print.gap = 2,	quote = FALSE)
cat("Alpha (exp):\n")
print.default(format(alpha, digits = digits), print.gap = 2, quote = FALSE)
cat("\n")
cat("Fixed effects:\n")
print.default(format(B, digits = digits), print.gap = 2, quote = FALSE)
cat("\n")
cat("Asymmetry:\n")
print.default(format(mu, digits = digits), print.gap = 2, quote = FALSE)
cat("\n")
cat("Scale matrix of the error:\n")
print.default(format(as.matrix(Sigma), digits = digits), quote = FALSE)
cat("\n")
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat("\n")
invisible(x)
}

update_par <- function(x, conditional){

K <- x$K

# no change
if(x$conditional == conditional) return(x$par)

# conditional to marginal
if(x$conditional & !conditional){
	Sigma <- buildSigma(x$gamma, K = K, conditional = x$conditional)
	var_sin <- Sigma[K+2,K+2]
	Sigma <- Sigma/var_sin
	gamma <- extractSigma(Sigma, K = K, conditional = conditional)
	val <- c(x$b, gamma, x$muw, x$tau)
}

# marginal to conditional
if(!x$conditional & conditional){
	stop("Can't do it")
}

return(val)
}

# Get marginal densities from cglr object

dtheta.cglr <- function(object, data = NULL, r = NULL){

K <- object$K

if(is.null(data)){
	theta <- object$theta
	mmx <- object$mmx
	mmu <- object$mmu
	mma <- object$mma
} else {
	theta <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[3]]))), data)
	mmx <- model.matrix(as.formula(paste("~ ", deparse(object$formula[[3]]))), data)
	mmu <- model.matrix(object$mu.formula, data)
	mma <- model.matrix(object$alpha.formula, data)
}

n <- length(theta)
if(object$conditional & is.null(r)){
	r <- object$r
}

if(object$conditional){
	if(length(r) != n) warning("Length of r (", length(r), ") is not the same as other variables (", n, ")")
}

m <- mmx %*% coef(object)
m <- m[,-c(1:K),drop=FALSE]
sigma <- getVarCov(object, scale = TRUE)
sigmas <- sigma[-c(1:K),-c(1:K)]
tau <- object$tau

if(single_gq <- all(colnames(mma) == "(Intercept)")){
        alpha <- exp(tau[1])
        alpha <- min(c(alpha, 1e10))
} else {
        xtau <- mma %*% tau
        alpha <- exp(xtau)
        alpha <- pmin(alpha, 1e10)
}

if(object$conditional){
	#tmp <- try(dcpgl(theta[j], m = m[i,], sigma = sigmas, shape = alpha[i], r = r[i]), silent = TRUE)
} else {
	val <- try(C_dpgl_vv(theta = as.vector(theta), m = as.matrix(m), sigma = as.matrix(sigmas), alpha = as.vector(alpha), Hv = as.integer(object$control$Hv), logd = as.logical(FALSE)), silent = TRUE)
}
return(val)
}

dw.cglr <- function(object, variable, data = NULL){

K <- object$K

if(is.null(data)){
	w <- object$w
	mmx <- object$mmx
	mmu <- object$mmu
	mma <- object$mma
} else {
	w <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[2]]))), data)
	mmx <- model.matrix(as.formula(paste("~ ", deparse(object$formula[[3]]))), data)
	mmu <- model.matrix(object$mu.formula, data)
	mma <- model.matrix(object$alpha.formula, data)
}
n <- nrow(w)
if(ncol(w) != K) stop("w must have K =", K, "columns")

m <- mmx %*% coef(object)
m <- m[,c(1:K),drop=FALSE]
muw <- mmu %*% matrix(object$muw, ncol = K)
sigma <- getVarCov(object, scale = TRUE)
sigmaw <- sigma[1:K,1:K,drop=FALSE]
alpha <- exp(mma %*% object$tau)
alpha <- pmin(alpha, 1e10)

if(variable == "theta") stop("Use dtheta.cglr or dens.cglr(..., variable = 'theta')")
sel <- match(variable, object$response)
if(is.na(sel)) stop("variable invalid")

val <- rep(NA, n)
for(i in 1:n){
	z <- as.vector(w[i,] - m[i,])
	val[i] <- try(nlmm::dmgl(z[sel], mu = as.vector(muw[i,sel]), sigma = sigmaw[sel,sel,drop=FALSE], shape = alpha[i]), silent = TRUE)
}
return(val)
}

dens.cglr <- function(object, variable = "theta", data = NULL){

K <- object$K

g <- function(xnsel, xsel, sel, nsel, K, m, sigma, muw, alpha){
	x <- rep(0, K+1)
	x[sel] <- xsel
	x[nsel] <- xnsel
	w <- x[1:K]
	theta <- x[K+1]
	ans <- C_loglik_proj_cglr(w = matrix(w, nrow = 1), theta = as.vector(theta), xb = matrix(m, nrow = 1), sigma = as.matrix(sigma), muw = matrix(muw, nrow = 1), alpha = as.vector(alpha), Hv = as.integer(20))
	exp(ans)
}

if(is.null(data)){
	w <- object$w
	theta <- object$theta
	mmx <- object$mmx
	mmu <- object$mmu
	mma <- object$mma
} else {
	w <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[2]]))), data)
	theta <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[3]]))), data)
	mmx <- model.matrix(as.formula(paste("~ ", deparse(object$formula[[3]]))), data)
	mmu <- model.matrix(object$mu.formula, data)
	mma <- model.matrix(object$alpha.formula, data)
}
n <- nrow(w)
if(ncol(w) != K) stop("w must have K =", K, "columns")

m <- mmx %*% coef(object)
muw <- mmu %*% matrix(object$muw, ncol = K)
sigma <- getVarCov(object, scale = TRUE)
alpha <- exp(mma %*% object$tau)
alpha <- pmin(alpha, 1e10)

sel <- match(variable, c(object$response, "theta"))
if(is.na(sel)) stop("variable invalid")
nsel <- setdiff(1:(K+1), sel)

y <- cbind(w, theta)
low <- apply(y, 2, min)
up <- apply(y, 2, max)
low[K+1] <- 0
up[K+1] <- 2*pi
low <- low[nsel]
up <- up[nsel]

val <- rep(NA, n)
for(i in 1:n){
	val[i] <- try(cubature::pcubature(f = g, lowerLimit = low, upperLimit = up, xsel = y[i,sel], sel = sel, nsel = nsel, K = K, m = m[i,], sigma = sigma, muw = muw[i,], alpha = alpha[i], tol = 1e-03, vectorInterface = FALSE)$integral, silent = TRUE)
}
return(val)
}

# Get joint densities from cglr object

dcglr <- function(object, data, log = FALSE){

K <- object$K
Hv <- object$control$Hv

if(is.null(data)){
	w <- object$w
	theta <- object$theta
	mmx <- object$mmx
	mmu <- object$mmu
	mma <- object$mma
} else {
	w <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[2]]))), data)
	theta <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[3]]))), data)
	mmx <- model.matrix(as.formula(paste("~ ", deparse(object$formula[[3]]))), data)
	mmu <- model.matrix(object$mu.formula, data)
	mma <- model.matrix(object$alpha.formula, data)
}

# get parameters
xb <- mmx %*% coef(object)
sigma <- buildSigma(object$gamma, K = K, conditional = object$conditional)
xmu <- mmu %*% matrix(object$muw, ncol = K)
tau <- coef(object, which = "tau")
if(all(colnames(mma) == "(Intercept)")){
	alpha <- exp(tau[1])
	alpha <- min(c(alpha, 1e10))
} else {
	xtau <- mma %*% tau
	alpha <- exp(xtau)
	alpha <- pmin(alpha, 1e10)
}

if(object$conditional){
	y <- cbind(w, r*cos(theta), r*sin(theta))
	val <- C_loglik_cond_cglr(y = as.matrix(y), xb = as.matrix(xb), sigma = as.matrix(sigma), mu = as.matrix(cbind(xmu, 0, 0)), alpha = as.vector(alpha), Hv = as.integer(Hv))
} else {
	val <- C_loglik_proj_cglr(w = as.matrix(w), theta = as.vector(theta), xb = as.matrix(xb), sigma = as.matrix(sigma), muw = as.matrix(xmu), alpha = as.vector(alpha), Hv = as.integer(Hv))
}

if(!log) val <- exp(val)

return(val)
}

# Get mode of joint density from cglr object

mode_cglr <- function(object, data = NULL){

	g <- function(x, xb, sigma, muw, alpha, Hv){
		K <- length(x) - 1
		val <- C_loglik_proj_cglr(w = matrix(x[1:K], nrow = 1), theta = x[K+1], xb = matrix(xb, nrow = 1), sigma = as.matrix(sigma), muw = matrix(muw, nrow = 1), alpha = as.vector(alpha), Hv = as.integer(Hv))
		return(-sum(val))
	}


if(is.null(data)){
	w <- object$w
	theta <- object$theta
	mmx <- object$mmx
	mmu <- object$mmu
	mma <- object$mma
} else {
	w <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[2]]))), data)
	theta <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[3]]))), data)
	mmx <- model.matrix(as.formula(paste("~ ", deparse(object$formula[[3]]))), data)
	mmu <- model.matrix(object$mu.formula, data)
	mma <- model.matrix(object$alpha.formula, data)
}


# dimensions
n <- length(theta)
K <- object$K
Hv <- object$control$Hv

# get parameters
xb <- mmx %*% coef(object)
sigma <- buildSigma(object$gamma, K = K, conditional = object$conditional)
xmu <- mmu %*% matrix(object$muw, ncol = K)
tau <- coef(object, which = "tau")
xtau <- mma %*% tau
alpha <- exp(xtau)
alpha <- pmin(alpha, 1e10)

val <- matrix(NA, n, K+1)
for(i in 1:n){
	tmp <- try(optim(par = c(w[i,],theta[i]), fn = g, xb = xb[i,], sigma = sigma, muw = xmu[i,], alpha = alpha[i], Hv = Hv), silent = TRUE)
	if(!inherits(tmp, "try-error")) val[i,] <- tmp$par
}

val
}

mode_cglr_theta <- function(object, data = NULL){

	g <- function(x, m, sigma, alpha, Hv){
		k <- length(x)
		val <- rep(NA, k)
		for(j in 1:k){
			val[j] <- C_dpgl(theta = as.numeric(x[j]), m = as.vector(m), sigma = as.matrix(sigma), alpha = as.numeric(alpha), Hv = as.integer(Hv), logd = as.logical(TRUE))
		}
		return(-val)
	}
	


if(is.null(data)){
	theta <- object$theta
	mmx <- object$mmx
	mmu <- object$mmu
	mma <- object$mma
} else {
	theta <- model.matrix(as.formula(paste("~ 0 +", deparse(object$formula[[2]][[3]]))), data)
	mmx <- model.matrix(as.formula(paste("~ ", deparse(object$formula[[3]]))), data)
	mmu <- model.matrix(object$mu.formula, data)
	mma <- model.matrix(object$alpha.formula, data)
}

# dimensions
n <- length(theta)
K <- object$K
Hv <- object$control$Hv

# get parameters
m <- mmx %*% coef(object)
sigma <- buildSigma(object$gamma, K = K, conditional = object$conditional)
tau <- coef(object, which = "tau")
xtau <- mma %*% tau
alpha <- exp(xtau)
alpha <- pmin(alpha, 1e10)

val <- rep(NA, n)
for(i in 1:n){
	tmp <- try(optimize(f = g, interval = c(0,2*pi), m = m[i,(K+1):(K+2)], sigma = sigmas, alpha = alpha[i], Hv = Hv), silent = TRUE)
	if(!inherits(tmp, "try-error")) val[i] <- tmp$minimum
}

val
}

#########################################################
### Fit PN and PGL distributions
#########################################################

# Loglikelihood for PN and PGL densities
loglik_pn <- function(eta, y){ 

if(any(is.nan(eta))) return(Inf)

m <- eta[1:2]
tau <- exp(eta[3])
rho <- tanh(eta[4])
sigma <- matrix(c(tau^2, rho*tau, rho*tau, 1), 2, 2)

val <- -sum(dpn(y, m = m, sigma = sigma, log = TRUE))
return(val)
}

loglik_pgl <- function(eta, y, Hv){

if(any(is.nan(eta))) return(Inf)

m <- eta[1:2]
tau <- exp(eta[3])
rho <- tanh(eta[4])
sigma <- matrix(c(tau^2, rho*tau, rho*tau, 1), 2, 2)
alpha <- exp(eta[5])

val <- dpglC(as.vector(y), m = as.vector(m), sigma = as.matrix(sigma), shape = alpha, Hv = Hv, log = TRUE)

return(-sum(val))
}

# Fit PN and PGL densities
pn.fit <- function(y, optimizer = "optim", method = "Nelder-Mead", maxit = 500, verbose = FALSE, keep.data = FALSE){

eta <- rep(0, 4)

if(optimizer == "nlminb"){
        fit <- do.call(nlminb, args = list(objective = loglik_pn, start = eta, y = y, control = list(iter.max = maxit, trace = as.numeric(verbose))))
        cat(fit$message, "\n")
        fit$loglik <- -fit$objective
}

if(optimizer == "optim"){
        fit <- do.call(optim, args = list(fn = loglik_pn, par = eta, y = y, control = list(maxit = maxit, trace = as.numeric(verbose))))
        if(fit$convergence == 0) cat("Convergence successful", "\n")
        if(fit$convergence == 1) cat("Reached maximum number of iterations without convergence", "\n")
        fit$loglik <- -fit$value
}

fit$m <- fit$par[1:2]
tau <- exp(fit$par[3])
rho <- tanh(fit$par[4])
fit$sigma <- matrix(c(tau^2, rho*tau, rho*tau, 1), 2, 2)
if(keep.data) fit$y <- y
fit$nobs <- length(y)

class(fit) <- "pnfit"

return(fit)
}

pgl.fit <- function(y, optimizer = "optim", method = "Nelder-Mead", Hv = 20, maxit = 500, verbose = FALSE, keep.data = FALSE){

eta <- rep(0, 5)


if(optimizer == "nlminb"){
        fit <- do.call(nlminb, args = list(objective = loglik_pgl, start = eta, y = y, Hv = Hv, control = list(iter.max = maxit, trace = as.numeric(verbose))))
        cat(fit$message, "\n")
        fit$loglik <- -fit$objective
}

if(optimizer == "optim"){
        fit <- do.call(optim, args = list(fn = loglik_pgl, par = eta, y = y, Hv = Hv, control = list(maxit = maxit, trace = as.numeric(verbose))))
        if(fit$convergence == 0) cat("Convergence successful", "\n")
        if(fit$convergence == 1) cat("Reached maximum number of iterations without convergence", "\n")
        fit$loglik <- -fit$value
}

fit$m <- fit$par[1:2]
tau <- exp(fit$par[3])
rho <- tanh(fit$par[4])
fit$sigma <- matrix(c(tau^2, rho*tau, rho*tau, 1), 2, 2)
fit$alpha <- exp(fit$par[5])
if(keep.data) fit$y <- y
fit$nobs <- length(y)

class(fit) <- "pglfit"

return(fit)
}

print.pnfit <- function(x, digits = max(3, getOption("digits") - 3), ...){

cat("Projected normal distribution", "\n")
cat("Log-likelihood:\n")
print.default(format(x$loglik, digits = digits), print.gap = 2,	quote = FALSE)
cat("Location (theta):\n")
print.default(format(x$m, digits = digits), print.gap = 2, quote = FALSE)
cat("\n")
cat("Scale (sigma):\n")
print.default(format(as.matrix(x$sigma), digits = digits), quote = FALSE)
cat("\n")
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat("\n")
invisible(x)

}

print.pglfit <- function(x, digits = max(3, getOption("digits") - 3), ...){

cat("Projected generalized Laplace distribution", "\n")
cat("Log-likelihood:\n")
print.default(format(x$loglik, digits = digits), print.gap = 2,	quote = FALSE)
cat("Location (theta):\n")
print.default(format(x$m, digits = digits), print.gap = 2, quote = FALSE)
cat("\n")
cat("Scale (sigma):\n")
print.default(format(as.matrix(x$sigma), digits = digits), quote = FALSE)
cat("\n")
cat("Tails (alpha):\n")
print.default(format(x$alpha, digits = digits), print.gap = 2, quote = FALSE)
cat("\n")
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat("\n")
invisible(x)

}

# Get loglikelihood from pnfit and pglfit objects
logLik.pnfit <- function(object, ...){
val <- object$loglik
attr(val, "df") <- length(object$par)
attr(val, "rdf") <- object$nobs - length(object$par)
return(val)
}

logLik.pglfit <- function(object, ...){
val <- object$loglik
attr(val, "df") <- length(object$par)
attr(val, "rdf") <- object$nobs - length(object$par)
return(val)
}

# Get density from pnfit and pglfit objects
dens.pnfit <- function(object, y = NULL){

if(is.null(y)){
	if(!is.null(object$y)){
		y <- object$y
	} else stop("Provide values for 'y'")
}

val <- dpn(y, m = object$m, sigma = object$sigma)

list(x = y, y = val)

}

dens.pglfit <- function(object, y = NULL){

if(is.null(y)){
	if(!is.null(object$y)){
		y <- object$y
	} else stop("Provide values for 'y'")
}


val <- dpglC(y, m = object$m, sigma = object$sigma, shape = object$alpha, Hv = 40)

list(x = y, y = val)

}

# Get summary stats from pnfit and pglfit objects
stats.pnfit <- function(object){

C <- integrate(function(theta, m, sigma) cos(theta)*dpn(theta, m, sigma), lower = 0, upper = 2*pi, m = object$m, sigma = object$sigma)$val
S <- integrate(function(theta, m, sigma) sin(theta)*dpn(theta, m, sigma), lower = 0, upper = 2*pi, m = object$m, sigma = object$sigma)$val
thetabar <- atan2(S, C)
rl <- sqrt(sum(c(C, S)^2))

list(circ_mean = thetabar, result_length = rl, circ_sd = sqrt(-2*log(rl)))
}

stats.pglfit <- function(object){

C <- integrate(function(theta, m, sigma, shape) cos(theta)*dpglC(theta, m, sigma, shape), lower = 0, upper = 2*pi, m = object$m, sigma = object$sigma, shape = object$alpha)$val
S <- integrate(function(theta, m, sigma, shape) sin(theta)*dpglC(theta, m, sigma, shape), lower = 0, upper = 2*pi, m = object$m, sigma = object$sigma, shape = object$alpha)$val
thetabar <- atan2(S, C)
rl <- sqrt(sum(c(C, S)^2))

list(circ_mean = thetabar, result_length = rl, circ_sd = sqrt(-2*log(rl)))
}

#########################################################
### Fit GL distribution
#########################################################

# Fit GL density
loglik_gl <- function(eta, y, Hv){ 

if(any(is.nan(eta))) return(Inf)

if(!is.matrix(y)) stop("y must be a matrix n x d")
d <- ncol(y)
n <- nrow(y)
g <- d*(d+1)/2

theta <- eta[1:d]
gamma <- eta[(d+1):(d+g)]
sigma <- as.matrix(nlme::pdSymm(value = gamma))
mu <- eta[(d+g+1):(d+g+d)]
alpha <- exp(eta[d+g+d+1])

val <- C_dmgl_vv(y = as.matrix(y), theta = as.vector(theta), sigma = as.matrix(sigma), mu = as.vector(mu), alpha = as.numeric(alpha), as.integer(n), Hv = as.integer(Hv), logd = as.logical(TRUE))

return(-sum(val))
}

gl.fit <- function(y, optimizer = "optim", method = "Nelder-Mead", Hv = 20, maxit = 500, verbose = FALSE, keep.data = FALSE){

if(is.matrix(y)){
	d <- ncol(y)
	n <- nrow(y)
}
if(is.vector(y)){
	d <- 1
	n <- length(y)
	y <- matrix(y, ncol = 1)
}

eta <- rep(0, d*(d+5)/2 + 1)

if(optimizer == "nlminb"){
        fit <- do.call(nlminb, args = list(objective = loglik_gl, start = eta, y = y, Hv = Hv, control = list(iter.max = maxit, trace = as.numeric(verbose))))
        cat(fit$message, "\n")
        fit$loglik <- -fit$objective
}

if(optimizer == "optim"){
        fit <- do.call(optim, args = list(fn = loglik_gl, par = eta, y = y, Hv = Hv, control = list(maxit = maxit, trace = as.numeric(verbose))))
        if(fit$convergence == 0) cat("Convergence successful", "\n")
        if(fit$convergence == 1) cat("Reached maximum number of iterations without convergence", "\n")
        fit$loglik <- -fit$value
}

g <- d*(d+1)/2

gamma <- fit$par[(d+1):(d+g)]
fit$theta <- as.numeric(fit$par[1:d])
fit$sigma <- as.matrix(nlme::pdSymm(value = gamma))
if(d == 1) fit$sigma <- sqrt(as.numeric(fit$sigma))
fit$mu <- as.numeric(fit$par[(d+g+1):(d+g+d)])
fit$alpha <- exp(fit$par[d+g+d+1])
if(keep.data) fit$y <- y
fit$nobs <- n
fit$dim <- d

class(fit) <- "glfit"

return(fit)
}

print.glfit <- function(x, digits = max(3, getOption("digits") - 3), ...){

cat("Generalized Laplace distribution", "\n")
cat("Log-likelihood:\n")
print.default(format(x$loglik, digits = digits), print.gap = 2,	quote = FALSE)
cat("Location (theta):\n")
print.default(format(x$theta, digits = digits), print.gap = 2, quote = FALSE)
cat("\n")
cat("Scale (sigma):\n")
print.default(format(as.matrix(x$sigma), digits = digits), quote = FALSE)
cat("\n")
cat("Asymmetry (mu):\n")
print.default(format(x$mu, digits = digits), print.gap = 2, quote = FALSE)
cat("\n")
cat("Tails (alpha):\n")
print.default(format(x$alpha, digits = digits), print.gap = 2, quote = FALSE)
cat("\n")
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat("\n")
invisible(x)

}

# Get GL density from glfit objects
dens.glfit <- function(object, y = NULL){

if(is.null(y)){
	if(!is.null(object$y)){
		y <- object$y
	} else stop("Provide values for 'y'")
}

theta <- object$theta
sigma <- object$sigma
if(object$dim == 1) sigma <- sigma^2
mu <- object$mu
alpha <- object$alpha
n <- length(y)

val <- C_dmgl_vv(y = as.matrix(y), theta = as.vector(theta), sigma = as.matrix(sigma), mu = as.vector(mu), alpha = as.numeric(alpha), as.integer(n), Hv = as.integer(40), logd = as.logical(FALSE))

list(x = y, y = val)

}

# Get mean and var from glfit objects
mean.glfit <- function(object, ...){
object$theta + object$mu*object$alpha
}

vcov.glfit <- function(object, ...){
if(object$dim == 1) object$sigma <- object$sigma^2
(object$sigma + tcrossprod(object$mu))*object$alpha
}

#########################################################
### Projected normal
#########################################################

# Joint density of polar coordinates of multivariate normal
dmvnormpol <- function(theta, r, m, sigma, log = FALSE){

# theta is the angle of the response - can be a matrix (k - 1) x n
# r is the length - a vector n x 1
# m is the location parameter - must be a vector of length k
# sigma is the k x k scale parameter of the multivariate normal in rectangular coordinates

flag <- inherits(theta, "array")

if(flag){
	k <- nrow(theta) + 1
	n <- ncol(theta)
} else {
	k <- length(theta) + 1
	n <- 1
	theta <- matrix(theta, k - 1, n)
}

if(length(r) != n) stop("r must be of length ", n) 
if(length(m) != k) stop("m must be of length ", k)
if(!all(dim(sigma) == k)) stop("dimension of sigma must be ", k, " by ", k)

z <- t(SphericalCubature::polar2rect(r = r, phi = theta)) # matrix n x k
e <- sweep(z, 2, m, "-")

ans <- mvtnorm::dmvnorm(e, mean = rep(0, k), sigma = sigma, log = TRUE) + (k-1)*log(r)

if(!log) ans <- exp(ans)

return(ans)

}

# vector-valued dmvnormpol
dmvnormpol_vv <- function(x, r, m, sigma){
		matrix(apply(x, 2, function(z, r, m, s) dmvnormpol(theta = z, r = r, m = m, sigma = sigma), r = r, m = m, s = sigma), ncol = ncol(x))
}

# Marginal density of angle from multivariate normal in polar coordinates (projected normal)
dpn <- function(theta, m, sigma, log = FALSE){

# theta is the angle of the response - can be a matrix (k - 1) x n
# m is the location parameter - must be a vector of length k
# sigma is the k x k scale parameter of the multivariate normal in rectangular coordinates

flag <- inherits(theta, "array")

if(flag){
	k <- nrow(theta) + 1
	n <- ncol(theta)
} else {
	k <- 2
	n <- length(theta)
	theta <- matrix(theta, k - 1, n)
}

if(!is.matrix(sigma)) stop("sigma must be a matrix")
if(!all(dim(sigma) == k)) stop("dimension of sigma must be ", k, " by ", k)
if(length(m) != k) stop("m must be of length ", k)

z <- t(SphericalCubature::polar2rect(r = rep(1, n), theta)) # matrix n x k
m <- matrix(m, k, 1) # matrix k x 1

ans <- rep(NA, n)
for(i in 1:n){
	ans[i] <- dmvpn_arma(as.vector(z[i,]), as.vector(m), as.matrix(sigma), as.logical(TRUE))
}

if(!log) ans <- exp(ans)
return(ans)
}

# Marginal density of length from multivariate normal in polar coordinates (Weil) using numerical cubature
dnweil <- function(r, m, sigma, log = FALSE){

n <- length(r)

if(!is.matrix(sigma)) stop("sigma must be a matrix")
if(dim(sigma)[1] != dim(sigma)[2]) stop("sigma must be square")
k <- dim(sigma)[1]
if(length(m) != k) stop("m must be a vector of length", k)

#last angle is between 0 and 2pi
ll <- c(rep(-pi/2, k - 2), -pi)
ul <- c(rep(pi/2, k - 2), pi)

# vectorized function is about the same speed as non-vectorized
# pcubature is faster than hcubature up to k = 3

ans <- rep(NA, n)

if(k < 4){
	for(i in 1:n){
		ans[i] <- cubature::pcubature(f = dmvnormpol_vv, lowerLimit = ll, upperLimit = ul, tol = 1e-3, vectorInterface = TRUE, r = r[i], m = m, sigma = sigma)$integral
		#ans <- cubature::pcubature(f = dmvnormpol, lowerLimit = ll, upperLimit = ul, tol = 1e-3, vectorInterface = FALSE, r = r[i], m = m, sigma = sigma)$integral
	}
} else {
	for(i in 1:n){
		ans[i] <- cubature::hcubature(f = dmvnormpol_vv, lowerLimit = ll, upperLimit = ul, tol = 1e-3, vectorInterface = TRUE, r = r[i], m = m, sigma = sigma)$integral
	}
}

if(log) ans <- log(ans)
return(ans)
}

# Conditional density of angle from multivariate normal in polar coordinates using numerical cubature
dcpn <- function(theta, m, sigma, r = 1, log = FALSE){

if(!is.matrix(sigma)) stop("sigma must be a matrix")
if(dim(sigma)[1] != dim(sigma)[2]) stop("sigma must be square")
k <- dim(sigma)[1]
if(length(m) != k) stop("m must be of length ", k)

theta <- matrix(theta, k - 1)
n <- ncol(theta)

if(length(r) != n){
	warning("theta (", n, ") and r (", length(r), ") have different dimensions. Recycling first value of r")
	r <- rep(r[1], n)
}

Cr <- dnweil(r = r, m = m, sigma = sigma)
ans <- dmvnormpol(theta = theta, r = r, m = m, sigma = sigma)/Cr

if(log) ans <- log(ans)

ans
}

# Conditional density of length from multivariate normal in polar coordinates
dcnweil <- function(r, m, sigma, theta = 0, log = FALSE){

if(!is.matrix(sigma)) stop("sigma must be a matrix")
if(dim(sigma)[1] != dim(sigma)[2]) stop("sigma must be square")
k <- dim(sigma)[1]
if(length(m) != k) stop("m must be of length ", k)

theta <- matrix(theta, k - 1)
n <- length(r)

if(length(r) != n){
	warning("theta (", ncol(theta), ") and r (", n, ") have different dimensions. Recycling first column of theta")
	r <- matrix(theta[,1], ncol = n)
}

Ctheta <- dpn(theta = theta, m = m, sigma = sigma, log = FALSE)
ans <- dmvnormpol(theta = theta, r = r, m = m, sigma = sigma, log = FALSE)/Ctheta

if(log) ans <- log(ans)

ans
}

# Summary stats for PN
stats.pn <- function(m, sigma){

C <- integrate(function(theta, m, sigma) cos(theta)*dpn(theta, m, sigma), lower = 0, upper = 2*pi, m = m, sigma = sigma)$val
S <- integrate(function(theta, m, sigma) sin(theta)*dpn(theta, m, sigma), lower = 0, upper = 2*pi, m = m, sigma = sigma)$val
thetabar <- atan2(S, C)
rl <- sqrt(sum(c(C, S)^2))

list(circ_mean = thetabar, result_length = rl, circ_sd = sqrt(-2*log(rl)))

}

#########################################################
### Projected generalized Laplace
#########################################################

# Joint density of polar coordinates of multivariate GL
dmglpol <- function(theta, r, m, sigma, shape, log = FALSE){

n <- length(theta)

if(length(r) != n) stop("theta and r must have the same length") 
if(length(m) != 2) stop("length of m must be 2")
if(!all(dim(sigma) == 2)) stop("dimension of sigma must be 2 by 2")

ans <- rep(0, n)
for(i in 1:n){
	z <- t(SphericalCubature::polar2rect(r = r[i], phi = theta[i])) # matrix 2 x 1
	e <- matrix(z - m)
	ans[i] <- as.numeric(nlmm::dmgl(e, mu = rep(0, 2), sigma = sigma, shape = shape, log = TRUE) + log(r[i])) # multiply by Jacobian r
}

if(!log) ans <- exp(ans)

return(ans)

}

# Marginal density of angle from multivariate GL in polar coordinates (projected GL)
dpgl <- function(theta, m, sigma, shape, log = FALSE){

n <- length(theta)

if(!is.matrix(sigma)) stop("sigma must be a matrix")
if(!all(dim(sigma) == 2)) stop("dimension of sigma must be 2 by 2")
if(length(m) != 2) stop("length of m must be 2")

# vectorize dmglpol wrt r
g <- function(r, theta, m, sigma, shape){
	val <- mapply(dmglpol, r = r, MoreArgs = list(theta = theta, m = m, sigma = sigma, shape = shape, log = FALSE))
	matrix(val, nrow = 1)
}

eps <- .Machine$double.eps

ans <- rep(0, n)
for(i in 1:length(theta)){
	# integrate slower than pcubature
	#ans[i] <- integrate(g, lower = eps, upper = Inf, theta = theta[i], m = m, sigma = sigma, shape = shape)$val
	ans[i] <- cubature::pcubature(f = g, lowerLimit = eps, upperLimit = Inf, tol = 1e-3, vectorInterface = TRUE, theta = theta[i], m = m, sigma = sigma, shape = shape)$integral
}

if(log) ans <- log(ans)

return(ans)

}

# Marginal density of angle from scale-mixture representation of multivariate GL (projected GL)
dpglC <- function(theta, m, sigma, shape, Hv = 20, log = FALSE){

n <- length(theta)

if(!is.matrix(sigma)) stop("sigma must be a matrix")
if(!all(dim(sigma) == 2)) stop("dimension of sigma must be 2 by 2")
if(length(m) != 2) stop("length of m must be 2")

ans <- rep(0, n)
for(i in 1:length(theta)){
	ans[i] <- C_dpgl(theta = as.numeric(theta[i]), m = as.vector(m), sigma = as.matrix(sigma), alpha = as.numeric(shape), Hv = as.integer(Hv), logd = as.logical(log))
}

return(ans)

}

# CDF of theta
ppgl <- function(theta, m, sigma, shape, log = FALSE){

n <- length(theta)

if(!is.matrix(sigma)) stop("sigma must be a matrix")
if(!all(dim(sigma) == 2)) stop("dimension of sigma must be 2 by 2")
if(length(m) != 2) stop("length of m must be 2")

eps <- .Machine$double.eps

ans <- rep(0, n)
for(i in 1:n){
	ans[i] <- integrate(dpgl, lower = eps, upper = theta[i], m = m, sigma = sigma, shape = shape, log = FALSE)$value
}

if(log) ans <- log(ans)

return(ans)

}

# QF of theta
qpgl <- function(p, m, sigma, shape){

n <- length(p)

if(!is.matrix(sigma)) stop("sigma must be a matrix")
if(!all(dim(sigma) == 2)) stop("dimension of sigma must be 2 by 2")
if(length(m) != 2) stop("length of m must be 2")

g <- function(x, p, m, sigma, shape){
	(p - ppgl(x, m, sigma, shape, log = FALSE))^2
}

eps <- .Machine$double.eps

ans <- rep(0, n)
for(i in 1:n){
	ans[i] <- optimize(f = g, lower = 0, upper = 2*pi, p = p[i], m = m, sigma = sigma, shape = shape)$minimum
}

return(ans)

}

# Marginal density of length from multivariate GL in polar coordinates (Weil) using numerical cubature
dglweil <- function(r, m, sigma, shape, log = FALSE){

n <- length(r)

if(!is.matrix(sigma)) stop("sigma must be a matrix")
if(!all(dim(sigma) == 2)) stop("dimension of sigma must be 2 by 2")
if(length(m) != 2) stop("length of m must be 2")

# vectorize dmglpol wrt theta
g <- function(theta, r, m, sigma, shape){
	val <- mapply(dmglpol, theta = theta, MoreArgs = list(r = r, m = m, sigma = sigma, shape = shape, log = FALSE))
	matrix(val, nrow = 1)
}

eps <- .Machine$double.eps

ans <- rep(0, n)
for(i in 1:n){
	#ans[i] <- integrate(g, lower = eps, upper = 2*pi, r = r[i], m = m, sigma = sigma, shape = shape)$val
	ans[i] <- pcubature(f = g, lowerLimit = eps, upperLimit = 2*pi, tol = 1e-3, vectorInterface = TRUE, r = r[i], m = m, sigma = sigma, shape = shape)$integral
}

if(log) ans <- log(ans)

return(ans)

}

# CDF
pglweil <- function(r, m, sigma, shape, log = FALSE){

n <- length(r)

if(!is.matrix(sigma)) stop("sigma must be a matrix")
if(!all(dim(sigma) == 2)) stop("dimension of sigma must be 2 by 2")
if(length(m) != 2) stop("length of m must be 2")

# vectorize dglweil wrt r
g <- function(r, m, sigma, shape){
	val <- mapply(dglweil, r = r, MoreArgs = list(m = m, sigma = sigma, shape = shape, log = FALSE))
	matrix(val, nrow = 1)
}

eps <- .Machine$double.eps
ans <- rep(0, n)
for(i in 1:n){
	#ans[i] <- integrate(dglweil, lower = eps, upper = r, m = m, sigma = sigma, shape = shape, log = FALSE)$value
	ans[i] <- pcubature(f = g, lowerLimit = eps, upperLimit = r[i], tol = 1e-3, vectorInterface = TRUE, m = m, sigma = sigma, shape = shape)$integral
}

if(log) ans <- log(ans)

return(ans)

}

# PDF of r analytic (equation 5)
dglweil2 <- function(x, shape, d = 2, log = FALSE){

if(any(x <= 0)) stop("x must be strictly positive")

N <- 2*x^(d/2+shape-1)*besselK(x=sqrt(2)*x, nu = d/2-shape)
D <- sqrt(2)^(d/2-2+shape)*gamma(shape)*gamma(d/2)

ans <- N/D

if(log) ans <- log(ans)

return(ans)

}

# Conditional density of angle from multivariate GL in polar coordinates using numerical cubature
dcpgl <- function(theta, m, sigma, shape, r = 1, log = FALSE){

if(length(r) != length(theta)){
	warning("theta (", length(theta), ") and r (", length(r), ") have different lengths. Recycling first value of r")
	r <- rep(r[1], length(theta))
}

if(all(m == 0)){
	Cr <- dglweil2(x = r, s = shape)
} else {
	Cr <- dglweil(r = r, m = m, sigma = sigma, shape = shape)
}

ans <- dmglpol(theta = theta, r = r, m = m, sigma = sigma, shape = shape)/Cr

if(log) ans <- log(ans)

ans
}

# Conditional density of length from multivariate GL in polar coordinates
dcglweil <- function(r, m, sigma, shape, theta = 0, log = FALSE){

if(length(r) != length(theta)){
	warning("r (", length(r), ") and theta (", length(theta), ") have different dimensions. Recycling first element of theta")
	theta <- matrix(theta[1], length(r))
}


Ctheta <- dpglC(theta = theta, m = m, sigma = sigma, shape = shape, log = FALSE)
ans <- dmglpol(theta = theta, r = r, m = m, sigma = sigma, shape = shape)/Ctheta

if(log) ans <- log(ans)

ans
}

# Summary stats for PGL
stats.pgl <- function(m, sigma, shape){

C <- integrate(function(theta, m, sigma, shape) cos(theta)*dpglC(theta, m, sigma, shape), lower = 0, upper = 2*pi, m = m, sigma = sigma, shape = shape)$val
S <- integrate(function(theta, m, sigma, shape) sin(theta)*dpglC(theta, m, sigma, shape), lower = 0, upper = 2*pi, m = m, sigma = sigma, shape = shape)$val
thetabar <- atan2(S, C)
rl <- sqrt(sum(c(C, S)^2))

list(circ_mean = thetabar, result_length = rl, circ_sd = sqrt(-2*log(rl)))

}
