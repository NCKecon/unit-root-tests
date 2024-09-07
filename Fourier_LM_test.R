
# /*
#   **> fourier_lm
# **  Reference:
#   **  Enders, W., and Lee, J. (2012),
# **  "A Unit Root Test Using a Fourier Series to Approximate Smooth Breaks"
# **  Oxford Bulletin of Economics and Statistics,74,4(2012),574-599.
# **
#   **  Format:  { LMk, k, p, cv[k, .] } = Fourier_LM(y[, pmax, kmax, ic]);
#   **
#     **  Input:   Y	     -  Nx1 matrix, data,
#   **
#     **           pmax    -  Optional, maximum number of lags for Ds; 0=no lags. Default = 8.
#     **
#       **           kmax    -  Optional, maximumum number of single Fourier frequency. Default = 5.
#       **                      (upper bound is 5)
#       **
#         **           ic      -  Optional, information criterion. Default = 3.:
#           **                      1=Akaike
#           **                      2=Schwarz
#           **                      3=t-stat significance
#           **
#             **  Output:  LMk     -  LM(k) statistic
#           **
#             **           k       -  Number of singlefrequency
#           **
#             **			 p       -  number of lags selected by chosen information criterion
#           **
#             **			 cv      -  1%, 5%, 10% critical values for the chosen model and k
#           **  Author: Saban Nazlioglu
#           -- The author makes no performance guarantees.
#           -- for public non-commercial use only.
#           -- for any bugs, please send e-mail to snazlioglu@pau.edu.tr
#           **
#             */
#             
#             /*03 February 2015 */
# This version of R code was written by Warren Kelly - original code - Saban Nazlioglu 


getFourierDeterministic <- function(y) {
  T <- length(y)  # Length of the series
  dy <- diff(y)  # First difference of y, length becomes T-1
  ly <- c(NA, y[-length(y)])  # Lagged y, same length as y
  dc <- rep(1, T)  # Constant vector
  dt <- 1:T  # Trend vector
  
  list(T = T, dy = dy, ly = ly, dc = dc, dt = dt)
}

runFourierOLS <- function(dep, z, k) {
  n <- nrow(z)
  m <- solve(t(z) %*% z)
  b <- m %*% t(z) %*% dep
  e <- dep - z %*% b
  ssr <- sum(e^2)
  sig2 <- ssr / (nrow(z) - ncol(z))
  se <- sqrt(diag(m) * sig2)
  # Calculate the t-statistic for the first coefficient
  taup <- b[1] / se[1]
  
  # Initialize variables
  aicp <- NA
  sicp <- NA
  tstatp <- NA
  if (k > 2) {  # Assuming Fourier terms are included
    # Calculate AIC and SIC
    aicp <- log(ssr / nrow(z)) + 2 * (k + 2) / nrow(z)
    sicp <- log(ssr / nrow(z)) + (ncol(z) + 2) * log(nrow(z)) / nrow(z)
    # t-statistic for the last coefficient
    tstatp <- b[length(b)] / se[length(b)]
  } else {
    # Log-likelihood
    LL <- -n / 2 * (1 + log(2 * pi) + log(ssr / n))
    # Calculate AIC
    aicp <- (2 * k - 2 * LL) / n
    # Calculate SIC
    sicp <- (k * log(n) - 2 * LL) / n
    # Absolute t-statistic for the last coefficient
    tstatp <- abs(b[length(b)] / se[length(b)])
  }
  
  # Return the computed values as a list
  list(taup = taup, aicp = aicp, sicp = sicp, tstatp = tstatp, ssrp = ssr)
}

getFourierDeterministic <- function(y) {
  T <- length(y)
  dy <- diff(y) 
  ly <- c(NA, y[-length(y)])
  dc <- rep(1, T)
  dt <- 1:T
  
  list(T = T, dy = dy, ly = ly, dc = dc, dt = dt)
}

detrendData <- function(y, z, dy, dsink, dcosk, p) {
  dz <- cbind(1, dsink, dcosk)
  b0 <- solve(t(dz) %*% dz) %*% t(dz) %*% dy
  psi_ <- y[1] - sum(dz[1, ] * b0)
  detrended_data <- y - psi_ - (z %*% b0)
  
  return(detrended_data)
}

detrendData <- function(y, z, dy, dsink, dcosk, p) {
  # Combine a column of ones, dsink, and dcosk into a matrix dz
  dz <- cbind(1, dsink, dcosk)
  
  # Estimate the coefficients b0 by regressing dy on dz
  if (p == 0) {
    b0 <- solve(t(dz) %*% dz) %*% t(dz) %*% dy
  } else {
    b0 <- solve(t(dz) %*% dz) %*% t(dz) %*% dy[-(1:(p+1))]
  }
  # Calculate psi_ (intercept adjustment)
  psi_ <- y[1] - sum(dz[1, ] * b0)
  # Return the detrended data: y minus psi_ and minus the product of z and b0
  detrended_data <- y - psi_ - (z %*% b0)
  
  return(detrended_data)
}

getFourierLMCrit <- function(T) {
  if (T <= 150) {
    return(matrix(c(-4.69, -4.10, -3.82,
                    -4.25, -3.57, -3.23,
                    -3.98, -3.31, -2.96,
                    -3.85, -3.18, -2.86,
                    -3.75, -3.11, -2.81), nrow = 5, byrow = TRUE))
  } else if (T <= 349) {
    return(matrix(c(-4.61, -4.07, -3.79,
                    -4.18, -3.55, -3.23,
                    -3.94, -3.30, -2.98,
                    -3.80, -3.18, -2.88,
                    -3.73, -3.12, -2.83), nrow = 5, byrow = TRUE))
  } else if (T <= 500) {
    return(matrix(c(-4.57, -4.05, -3.78,
                    -4.13, -3.54, -3.22,
                    -3.94, -3.31, -2.98,
                    -3.81, -3.19, -2.88,
                    -3.75, -3.14, -2.83), nrow = 5, byrow = TRUE))
  } else {
    return(matrix(c(-4.56, -4.03, -3.77,
                    -4.15, -3.54, -3.22,
                    -3.94, -3.30, -2.98,
                    -3.80, -3.19, -2.88,
                    -3.74, -3.13, -2.83), nrow = 5, byrow = TRUE))
  }
}

get_lag <- function(ic, pmax, aicp, sicp, tstatp) {
  # Initialize p
  p <- NA
  # Information Criterion: 1 = Akaike
  if (ic == 1) {
    p <- which.min(aicp)
  }
  # Information Criterion: 2 = Schwarz (BIC)
  if (ic == 2) {
    p <- which.min(sicp)
  }
  # Information Criterion: 3 = t-stat significance
  if (ic == 3) {
    j <- pmax + 1
    isw <- 0
    while (isw != 1) {
      if (abs(tstatp[j]) > 1.645 || j == 1) {
        p <- j
        isw <- 1
      }
      j <- j - 1
    }
    # Ensure p is valid
    if (is.na(p) || p == 0) {
      p <- 1
    }
  }
  return(p)
}

Fourier_LM <- function(y, pmax = 8, kmax = 5, ic = 3) {
  y <- as.vector(y)
  T <- length(y)
  
  # Initialize variables
  taup <- numeric(pmax + 1)
  aicp <- numeric(pmax + 1)
  sicp <- numeric(pmax + 1)
  tstatp <- numeric(pmax + 1)
  ssrp <- numeric(pmax + 1)
  ssrk <- numeric(kmax)
  tauk <- numeric(kmax)
  keep_p <- numeric(kmax)
  
  # Prepare data
  for (k in 1:kmax) {
    #k=1
    # Fourier Deterministic Setup
    fourier_det_y <- getFourierDeterministic(y)
    dt <- fourier_det_y$dt
    dy <- fourier_det_y$dy
    ly <- fourier_det_y$ly
    dc <- fourier_det_y$dc
    
    # Fourier Terms
    t <- 1:T
    sinp <- sin(2 * pi * k * t / T)
    cosp <- cos(2 * pi * k * t / T)
    dsink <- diff(sinp, 1)
    dcosk <- diff(cosp, 1)
    z <- cbind(dt, sinp, cosp)
    
    ylm <- detrendData(y, z, dy, dsink, dcosk, p = 0)
    
    new_det_trend <- getFourierDeterministic(ylm)
    #t <- new_det_trend$t
    dy <- new_det_trend$dy
    ly <- new_det_trend$ly
    dc <- new_det_trend$dc
    dt <- new_det_trend$dt
    
    # New Fourier Terms
    sink <- diff(sinp)
    cosk <- diff(cosp)
    
    for (p in 0:pmax) {
      
      # Trim `ldy`, `y1`, `sbt`, and `trnd` consistently
      # Create Regression Matrix z
      ldy <- if (p > 0) { embed(dy, p + 1)[, -1] } else {
        rep(0, (T-1))
      }
      
      sbt <- dc[-(1:(p + 1))]
      trnd <- dt[-(1:(p + 1))]
      
      if (p == 0) {
        y1 <- ly[-1]
        sbt <- dc[-1]
        dep <- dy
        sinp <- sink
        cosp <-cosk
        z <- cbind(y1, sbt, sinp, cosp)
      } else {
        y1 <- ly[-(1:(1+p))]
        sbt <- sbt
        dep <- dy[-(1:(p))]
        sinp <- sink[-(1:(p))]
        cosp <-cosk[-(1:(p))]
        z <- cbind(y1, sbt, sinp, cosp, ldy)
      }
      
      # OLS Regression
      F_OLS <- runFourierOLS(dep, z, k)
      ssrp[p + 1] <- F_OLS$ssr
      taup[p + 1] <- F_OLS$taup
      aicp[p + 1] <- F_OLS$aicp
      sicp[p + 1] <- F_OLS$sicp
      tstatp[p + 1] <- F_OLS$tstatp
    }
    
    # Select optimal lag
    p_opt <- get_lag(ic, pmax, aicp, sicp, tstatp)
    keep_p[k] <- p_opt
    ssrk[k] <- ssrp[p_opt]
    tauk[k] <- taup[p_opt]
  }
  
  # Find optimal frequency
  f <- which.min(ssrk)
  LMk <- tauk[f]
  opt_lag <- keep_p[f]
  
  # Critical values
  crit <- getFourierLMCrit(T)
  
  # Print results
  cat("Fourier LM Test (Enders & Lee, 2012)\n")
  cat("LM statistic:", LMk, "\n")
  cat("Optimal Frequency:", f, "\n")
  cat("Optimal Lag:", opt_lag, "\n")
  cat("Critical Values: 1%:", crit[f, 1], "5%:", crit[f, 2], "10%:", crit[f, 3], "\n")
  
  return(list(LMk = LMk, f = f, opt_lag = opt_lag, crit = crit[f, ]))
}


fourier_gls(y, model = 1, pmax = 4, fmax = 3, ic = 1)
