
library(GIGrvg)

# función de pérdida.
check_loss<-function(t,theta){
 ifelse(t<0,(theta-1)*t,theta*t)
}

# función de la log-apriori.
log_prior <- function(b, d, varB, tau, mu, pin, prior_dist, a_hyp, b_hyp, c_hyp, d_hyp) {
 
  log_p_b <- sum(dnorm(b, mean = 0, sd = sqrt(varB), log = TRUE))
  
  log_p_d <- sum(d * log(pin) + (1 - d) * log(1 - pin))
  
  if (length(b_hyp) != 1L) {
    stop("log_prior(): 'b_hyp' debe ser escalar (rate de la Gamma de tau). ",
         "Parece que se pasó el vector de coeficientes 'b'.")
  }
  
  if (prior_dist$class == "Inv-chi-square") {
    dfb <- prior_dist$dfb
    Sb <- prior_dist$Sb
    log_p_varB <- ((dfb/2) * log(Sb/2) - lgamma(dfb/2)) - (((dfb/2) + 1) * log(varB)) - (Sb / (2 * varB))
  } else if (prior_dist$class == "Half-Cauchy") {
    Sb <- prior_dist$Sb
    log_p_varB <- log(2) - log(pi) - log(Sb) - log(1 + (varB / Sb)^2)
  } else {
    stop("Distribución a priori no válida. Usar 'Inv-chi-square' o 'Half-Cauchy'.")
  }
  
  log_p_tau <- (a_hyp - 1) * log(tau) - b_hyp * tau
  
  log_p_mu <- 0
  
  log_p_pi <- dbeta(pin, shape1 = c_hyp, shape2 = d_hyp, log = TRUE)
  
  log_prior_total <- log_p_b+log_p_d+log_p_varB+log_p_tau+log_p_mu+log_p_pi
  
  return(log_prior_total)
}

# generar la función de la logLike aumentada.
log_likeA <- function(y, X, beta, mu, tau, v, e1, e2) 
{
  n <- length(y)
  v_i <- pmax(v, 1e-12)
  ytilde <- y - mu - as.vector(X %*% beta) - e1 * v_i
  
  ll <- (3*n/2) * log(tau) -
    n * log(e2) -
    0.5 * sum(log(v_i)) -
    (tau / (2 * e2^2)) * sum((ytilde^2) / v_i) -
    tau * sum(v_i)
  
  return(ll)
}

# función BRQRC.
BRQRC <- function(y, X, 
                  intercept = TRUE, 
                  theta = 0.50, 
                  prior_dist = list(class="Inv-chi-square", dfb = 5, Sb = NULL),
                  # prior_dist = list(class="Half-Cauchy", dfb = 1, Sb = 25),
                  a = 0.1, 
                  b = 0.1, 
                  probIn = 0.5,
                  counts = 2,   
                  R2 = 0.50,    
                  nIter = 1000,
                  burnIn = 500,
                  thin = 10,      
                  verbose = TRUE, 
                  save = FALSE,   
                  file_betas="betas.txt",
                  file_deltas="deltas.txt"
)      
{
  n <- length(y)
  p <- ncol(X)

  e1 <- (1 - 2*theta) / (theta*(1 - theta))
  e2 <- sqrt(2 / (theta*(1 - theta)))
  
  vy <- var(y)
  error <- y
  
  X2 <- X^2
  
  if(prior_dist$class%in%c("Inv-chi-square","Half-Cauchy"))
  {
    if(prior_dist$class=="Inv-chi-square")
    {
      
      message("Ha seleccionado la distribución Inv-chi-square\n")
      
      if(is.null(prior_dist$Sb))
      {
        x2 <- colSums(X2)
        sx <- colSums(X)
        sumMeanXSq <- sum((sx/n)^2)
        MSx <- sum(x2) / n-sumMeanXSq
        
        prior_dist$Sb<-R2*vy*(prior_dist$dfb+2)/MSx/probIn
        cat("Sb=",prior_dist$Sb,"\n")
      }
    }else{
      message("Ha seleccionado la distribución Half-Cauchy\n")
      
      if(is.null(prior_dist$dfb))
      {
        prior_dist$dfb = 1
        message("dfb se fijó a ", prior_dist$dfb,"\n")
      }
      if(is.null(prior_dist$Sb))
      {
        prior_dist$Sb = sd(y)/mean(apply(X,2,sd))
        message("Sb se fijó a ", prior_dist$Sb,"\n")
      }
    }
    
  }else{
    stop("class deber ser 'Inv-chi-square' o 'Half-Cauchy' \n")
  }
  
  varB <- 0.10            
  b_tau <- b
  
  b <- rep(0, p)
  d <- rep(1, p)
  D <- matrix(NA, nrow = nIter, ncol = p)
  beta <- b * d
  Beta <- matrix(NA, nrow = nIter, ncol = p)
  vVarB <- rep(NA, nIter)
  
  countsIn <- counts * probIn
  countsOut <- counts * (1 - probIn)
  
  PIN <- rep(NA, nIter)
  
  logLike <- rep(NA, nIter)
  logPrior <- rep(NA, nIter)
  logLikeA <- rep(NA, nIter)
  
  yHat <- matrix(NA, nrow = nIter, ncol = n)
  
  if(intercept)
  {
    Mu <- rep(NA, nIter)
    mu <- mean(y)
    error <- error - mu
  }else{
    mu <- 0
  }
  
  Tau <- rep(NA, nIter)
  tau <- 1                 
  
  v <- rexp(n = n)
  pin <- 0.5
  
  for(iter in 1:nIter)
  {
    time0 <- proc.time()[3]
    
    if(verbose)
    {
      cat("*********************************************************\n")
      cat("iter=",iter,"\n")
    }
    
    if(intercept)
    {
      ytilde <- error - e1*v + mu
      mean <- sum(ytilde/v)/sum(1/v)
      variance <- e2^2/(tau*sum(1/v))
      mu_new <- rnorm(1, mean = mean, sd = sqrt(variance))
      
      error <- error + (mu - mu_new)
      mu <- mu_new
      Mu[iter] <- mu
      
      if(verbose) cat("mu=",mu,"\n")
    }
    
    logPriorOdds <- log(pin / (1-pin))
    c1 <- (e2^2/tau) / varB
    
    for(j in 1:p)
    {	
      xj <- X[,j]
      ytilde <- error - e1*v + xj*beta[j]
      beta_old <- beta[j]
      
      tmp1 <- xj/v
      tmp2 <- tmp1*ytilde
      tmp3 <- tmp1*xj
      tmp4 <- sum(tmp2)
      
      logP <- tau/(2*e2^2) * (2*b[j] * tmp4 - b[j]^2*sum(tmp3)) + logPriorOdds
      
      u <- runif(1)
      logOdds <- log(u / (1-u))
          
      if(logP > logOdds)
      {
        d[j] <- 1
        g <- sum(X2[,j] / v) + c1
        b[j] <- rnorm(n = 1, mean = tmp4/g, sd = e2/sqrt(tau*g))
        beta[j] <- b[j]
      }else{    
        d[j] <- 0
        b[j] <- rnorm(n = 1, mean = 0, sd = sqrt(varB))
        beta[j] <- 0
      }
      
      error <- error + xj*(beta_old - beta[j])
      
    }
    Beta[iter,] <- beta
    D[iter,] <- d

    if(prior_dist$class=="Inv-chi-square"){
      varB <- (prior_dist$Sb + sum(b^2)) / rchisq(n = 1, df = prior_dist$dfb + p)
      vVarB[iter] <- varB
      
      if(verbose) cat("varB=",varB,"\n")
    }
    
    if (prior_dist$class=="Half-Cauchy") {
      
      W <- min(1 / rgamma(1, shape = 1, rate = (prior_dist$dfb/varB + 1/(prior_dist$Sb^2))), 1e3)
      
      varBtmp <- 1 / rgamma(1, shape = (p + prior_dist$dfb) / 2, rate = (sum(b^2) / 2) + (prior_dist$dfb / W))
      
      if(varBtmp < 1e-3 || varBtmp > 1e3) 
      {
        cat("warning: varB not updated\n")
        
        vVarB[iter]<-varB
        
      }else{
        
        vVarB[iter]<-varBtmp
        varB<-varBtmp
      }
      
      if(verbose) cat("varB=",varBtmp,"\n")
      
    }
    
    for(i in 1:n)
    {
      v[i] <- rgig(n=1,lambda=1/2,chi=tau*(error[i]/e2)^2,psi=tau*(e1/e2)^2+2*tau)
    }
    
    tmp5 <- b_tau+sum((error - e1*v)^2 / (2*e2^2*v) + v)
    tau <- rgamma(n = 1, shape = a+3*n/2, rate = tmp5)
    Tau[iter] <- tau
    
    if(verbose) cat("tau=",tau,"\n")
    
    sumIn <- sum(d)
    sumOut <- p - sumIn
    pin = rbeta(shape1 = sumIn + countsIn, shape2 = sumOut + countsOut, n = 1)
    PIN[iter] <- pin
    
    slog <- 0
    for(i in 1:n)
    {
      slog <- slog  + check_loss(t = error[i] * tau, theta = theta)
    }
    logLike[iter] <- (n * log(theta)) - (n * log(1-theta)) - (n * log(1 / tau)) - slog

    logPrior[iter] <- log_prior(b, d, varB, tau, mu, pin, prior_dist,
                                a, b_tau, countsIn, countsOut)
    
    logLikeA[iter] <- log_likeA(y, X, beta, mu, tau, v, e1, e2)
    
    
    yHat[iter, ] <- X %*% beta + if (intercept) mu else 0
    
    if(verbose) cat("pin=",pin,"\n")
    
    time1 <- proc.time()[3]
    
    if(verbose) cat("tiempo/iteración=", round(time1-time0,4)," segundos \n")
    
  }
  
  out = list()
  
  index <- seq(from = burnIn + thin, to = nIter, by = thin)
  
  if(save) write.table(Beta, file = file_betas, row.names=FALSE)
  if(save) write.table(D, file = file_deltas, row.names=FALSE)
  
  Beta <- Beta[index,]
  vVarB <- vVarB[index]
  
  out$beta <- colMeans(Beta)
  out$SD.beta <- as.vector(apply(Beta,2,sd))
  out$varB <- mean(vVarB)
  out$SD.varB <- sd(vVarB)
  
  if(intercept)
  {
    Mu = Mu[index]
    out$mu = mean(Mu)
    out$SD.mu = sd(Mu)
    error = error - out$mu
  }
  
  Tau <- Tau[index]
  out$tau <- mean(Tau)
  out$SD.tau <- sd(Tau)
  
  D <- D[index,]
  out$d <- colMeans(D)
  out$SD.d <- as.vector(apply(D,2,sd))
  PIN <- PIN[index]
  out$probIn <- mean(PIN)
  out$SD.probIn <- sd(PIN)
  out$numIter <- index
  out$logLike <- logLike
  out$logPrior <- logPrior
  out$logLikeA <- logLikeA
  out$logPP <- logLikeA + logPrior
  yHat <- yHat[index, ]
  out$yHat <- colMeans(yHat)
  out$SD.yHat <- as.vector(apply(yHat, 2, sd))

  return(out)
}
