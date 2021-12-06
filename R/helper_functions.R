##############################################################################
# Helper functions needed to run analyses in Sample Size Estimation using a 
# Latent Variable Model for Mixed Outcome Co-Primary, Multiple Primary and 
# Composite Endpoints
#
# Author: Martina McMenamin
#
#############################################################################



# Composite endpoint functions --------------------------------------------

#' @param mean.lat overall effect size on the composite using the latent variable model
#' @param mean.bin overall effect size on the composite using the standard binary method
#' @param var.lat estimated variance of effect using the latent variable model
#' @param var.bin variance of effect using the standard binary method
#' @param alpha significance level
#' @param beta 1-desired power 

samplesize <- function(mean.lat,
                       mean.bin,
                       var.lat,
                       var.bin,
                       alpha,
                       beta){
  
  n.lat=(var.lat*(qnorm(1-beta)-(qnorm(alpha)))^2)/(mean.lat)^2
  n.bin=(var.bin*(qnorm(1-beta)-(qnorm(alpha)))^2)/(mean.bin)^2
  
  return(c(n.lat,n.bin))
}


#' @param mean overall effect size on the composite 
#' @param var variance of effect (standard deviation squared)
#' @param alpha significance level
#' @param n overall sample size


powerfunc<-function(mean,var,alpha,n){
  
  power.lat=pnorm((mean/(sqrt(var/n)))-qnorm(1-alpha))
  
  return(power.lat)
}



# Family wise error rate for multiple primary endpoints -------------------

#' @param K number of endpoints
#' @param n sample size per arm
#' @param delta effect size for each outcome
#' @param Sigma covariance matrix for the outcomes
#' @param SD standard deviation of outcomes
#' @param rho correlation between endpoints 
#' @param adjust logical whether or not to implement Bonferroni correction
#' @param sig.level alpha level
#' @param power required power


FWER <- function (K, n = NULL, delta = NULL, Sigma, SD, rho, adjust, sig.level = 0.025, 
                  power = NULL, tol = .Machine$double.eps^0.25) 
{
  
  sig.level = ifelse(adjust == TRUE, sig.level/K, sig.level)
  
  if (is.null(power)) {
    std.effect <- delta/sqrt(diag(Sigma))
    z.alpha <- qnorm(1 - sig.level)
    crit.vals <- z.alpha - sqrt(n/2) * std.effect
    power <- 1 - pmvnorm(lower = -crit.vals, sigma = Sigma)
  }
  if (is.null(n)) {
    std.effect <- delta/sqrt(diag(Sigma))
    z.alpha <- qnorm(1 - sig.level)
    ssize.fct <- function(n, std.effect, z.alpha, Sigma, 
                          power) {
      crit.vals <- z.alpha - sqrt(n/2) * std.effect
      pmvnorm(lower = -crit.vals, sigma = Sigma) - 
        (1 - power)
    }
    n <- uniroot(ssize.fct, c(2, 1e+05), tol = tol, extendInt = "yes", 
                 std.effect = std.effect, z.alpha = z.alpha, Sigma.cor = Sigma, 
                 power = power)$root
  }
  NOTE <- "n is number in *each* group"
  METHOD <- "FWER for multiple primary endpoints"
  structure(list(n = n, delta = delta, SD = sqrt(diag(Sigma)), 
                 rho = Sigma[lower.tri(Sigma)], Sigma = Sigma, 
                 sig.level = sig.level, FWER = power, note = NOTE, method = METHOD), 
            class = "power.mpe.test")
}


