################################################################################
#### Zero-inflated Negative Binomial Models ####################################
################################################################################

## Hurdle Version ##############################################################

# log-likelihood function and derivatives

#' @rdname hurdNBlss
#' @export
hurdNBll <- function(y, eta1, eta2, eta3, deriv = 0) {

  ## ================================================================
  ## Hurdle zero-inflated negative binomial log likelihood
  ##
  ## Parameterisation:
  ##   mu  = exp(eta1)
  ##   p   = plogis(eta2)      -- EXACT P(Y=0), not a mixture weight
  ##   phi = exp(eta3)
  ##
  ##   P(Y=0)   = p
  ##   P(Y=y)   = (1-p) * f_NB(y|mu,phi) / (1 - f_NB(0|mu,phi))   (y>0)
  ##
  ## NB2:  Var(Y | NB) = mu + mu^2/phi
  ##
  ## deriv:
  ##   0 : likelihood only
  ##   1 : likelihood, first and second derivatives
  ##   2 : additionally third derivatives
  ##   3/4 : additionally fourth derivatives
  ##
  ## ----------------------------------------------------------------
  ## NOTE on structure vs ziNBll (the mixture-model sibling):
  ##
  ## The zero-count branch here is exact and trivial: ell_0 = log(p),
  ## depending on eta2 ONLY -- no posterior "q" weighting is needed, in
  ## contrast to the mixture model.
  ##
  ## The positive-count branch is
  ##   ell_+ = log(1-p) + A(y;eta1,eta3) - C(eta1,eta3)
  ## where A is the ordinary NB2 log-density (as in ziNBll) and
  ##   C = log(1 - f_NB(0|mu,phi))
  ## is the zero-truncation normalising constant. C is a composition
  ## C = h(A*) of h(x) = log(1-e^x) with the NB log-density evaluated
  ## at y=0 (A*, and its derivatives A1*..A3333* -- identical formulas
  ## to ziNBll's "starred" quantities). C's own derivatives follow from
  ## the (closed-form, tractable) derivatives of h via the multivariate
  ## Faa di Bruno formula; see the accompanying derivation.
  ##
  ## NOTE on El2 (expected Hessian):
  ## The (eta1,eta1), (eta1,eta3) blocks are FULLY closed-form (no
  ## intractable term appears in this hurdle model there); (eta1,eta2)
  ## and (eta2,eta3) are EXACTLY zero; (eta2,eta2) is the trivial
  ## -p(1-p). The (eta3,eta3) entry still requires the same intractable
  ## E_NB[psi_1(Y+phi)] term as ziNBll's Abar33 -- per current decision,
  ## this uses the OBSERVED curvature l2[,6] directly (as in ziNBll),
  ## not an expected-value approximation.
  ## ================================================================


  ## ------------------------------------------------
  ## Input checks
  ## ------------------------------------------------

  n <- length(y)

  if (length(eta1) != n ||
      length(eta2) != n ||
      length(eta3) != n) {
    stop("y, eta1, eta2 and eta3 must have equal lengths")
  }

  if (length(deriv) != 1L ||
      !is.finite(deriv) ||
      deriv != as.integer(deriv) ||
      deriv < 0 || deriv > 4) {
    stop("deriv must be an integer between 0 and 4")
  }

  if (any(!is.finite(y)) ||
      any(y < 0) ||
      any(y != floor(y))) {
    stop("y must contain finite non-negative integers")
  }

  deriv <- as.integer(deriv)


  ## ------------------------------------------------
  ## Stable inverse links
  ## ------------------------------------------------

  lo <- log(.Machine$double.xmin)
  hi <- log(.Machine$double.xmax)

  e1 <- pmin(pmax(eta1, lo), hi)
  e3 <- pmin(pmax(eta3, lo), hi)

  mu  <- exp(e1)
  phi <- exp(e3)

  p <- numeric(n)
  ip <- eta2 >= 0
  p[ip]  <- 1 / (1 + exp(-eta2[ip]))
  ep <- exp(eta2[!ip])
  p[!ip] <- ep / (1 + ep)

  om <- 1 - p


  ## ------------------------------------------------
  ## Allocate outputs
  ## ------------------------------------------------

  l <- numeric(n)

  l1 <- l2 <- l3 <- l4 <- El2 <- NULL

  if (deriv >= 1) {
    l1 <- matrix(0, n, 3L)
    l2 <- matrix(0, n, 6L)
    El2 <- matrix(0, n, 6L)
  }
  if (deriv > 1) l3 <- matrix(0, n, 10L)
  if (deriv > 3) l4 <- matrix(0, n, 15L)


  ## ------------------------------------------------
  ## Common NB quantities
  ## ------------------------------------------------

  D <- mu + phi
  log_phi_over_D <- log(phi) - log(D)

  Astar_g <- phi * log_phi_over_D
  w_g <- exp(Astar_g)

  ## ------------------------------------------------
  ## Row-wise starred (y=0) NB quantities and truncation-constant C
  ## derivatives. Needed for EVERY row (both the El2 block and the
  ## positive-count likelihood), since these depend only on mu, phi.
  ## ------------------------------------------------

  if (deriv >= 1) {

    a1_g  <- -mu * phi / D
    a3_g  <- phi * (log_phi_over_D + mu / D)
    a11_g <- -mu * phi^2 / D^2
    a13_g <- -mu^2 * phi / D^2
    a33_g <- a3_g + phi * mu^2 / D^2

    h1_g <- -w_g / (1 - w_g)
    h2_g <- -w_g / (1 - w_g)^2

    C1_g  <- h1_g * a1_g
    C3_g  <- h1_g * a3_g
    C11_g <- h2_g * a1_g^2 + h1_g * a11_g
    C13_g <- h2_g * a1_g * a3_g + h1_g * a13_g
    C33_g <- h2_g * a3_g^2 + h1_g * a33_g

    if (deriv > 1) {

      C4_g <- (mu - phi)^2 - 2 * mu * phi

      a111_g <- mu * phi^2 * (mu - phi) / D^3
      a113_g <- -mu * phi^2 * (mu - phi) / D^3 - mu * phi^2 / D^2
      a133_g <- -mu^2 * phi * (mu - phi) / D^3
      a333_g <- a3_g + 2 * phi * mu^3 / D^3

      h3_g <- -w_g * (1 + w_g) / (1 - w_g)^3

      C111_g <- h3_g * a1_g^3 + 3 * h2_g * a1_g * a11_g + h1_g * a111_g
      C113_g <- h3_g * a1_g^2 * a3_g +
        h2_g * (a11_g * a3_g + 2 * a1_g * a13_g) + h1_g * a113_g
      C133_g <- h3_g * a1_g * a3_g^2 +
        h2_g * (2 * a3_g * a13_g + a1_g * a33_g) + h1_g * a133_g
      C333_g <- h3_g * a3_g^3 + 3 * h2_g * a3_g * a33_g + h1_g * a333_g

      if (deriv > 3) {

        a1111_g <- -mu * phi^2 * C4_g / D^4
        a1113_g <- mu * phi^2 * C4_g / D^4 + mu * phi^2 * (mu - phi) / D^3
        a1133_g <- -mu * phi^2 * C4_g / D^4 - mu * phi^2 * (3 * mu - phi) / D^3
        a1333_g <- -mu^2 * phi * C4_g / D^4
        a3333_g <- a3_g + phi * mu^2 * (3 * mu^2 - 2 * mu * phi + phi^2) / D^4

        h4_g <- -w_g * (1 + 4 * w_g + w_g^2) / (1 - w_g)^4

        C1111_g <- h4_g * a1_g^4 + 6 * h3_g * a1_g^2 * a11_g +
          h2_g * (4 * a1_g * a111_g + 3 * a11_g^2) + h1_g * a1111_g
        C1113_g <- h4_g * a1_g^3 * a3_g +
          h3_g * (3 * a1_g * a3_g * a11_g + 3 * a1_g^2 * a13_g) +
          h2_g * (3 * a1_g * a113_g + a3_g * a111_g + 3 * a11_g * a13_g) +
          h1_g * a1113_g
        C1133_g <- h4_g * a1_g^2 * a3_g^2 +
          h3_g * (a1_g^2 * a33_g + 4 * a1_g * a3_g * a13_g + a3_g^2 * a11_g) +
          h2_g * (2 * a1_g * a133_g + 2 * a3_g * a113_g +
                    a11_g * a33_g + 2 * a13_g^2) +
          h1_g * a1133_g
        C1333_g <- h4_g * a1_g * a3_g^3 +
          h3_g * (3 * a1_g * a3_g * a33_g + 3 * a3_g^2 * a13_g) +
          h2_g * (3 * a3_g * a133_g + a1_g * a333_g + 3 * a33_g * a13_g) +
          h1_g * a1333_g
        C3333_g <- h4_g * a3_g^4 + 6 * h3_g * a3_g^2 * a33_g +
          h2_g * (4 * a3_g * a333_g + 3 * a33_g^2) + h1_g * a3333_g
      }
    }
  }


  ## ================================================================
  ## POSITIVE COUNTS: ell_+ = log(1-p) + A(y) - C
  ## ================================================================

  iz <- y == 0
  ipos <- !iz

  if (any(ipos)) {

    yy <- y[ipos]
    mm <- mu[ipos]
    pp <- p[ipos]
    oo <- om[ipos]
    ff <- phi[ipos]
    DD <- D[ipos]

    A <- lgamma(yy + ff) -
      lgamma(ff) -
      lgamma(yy + 1) +
      ff * log(ff / DD) +
      yy * log(mm / DD)

    Cval <- log(1 - w_g[ipos])

    l[ipos] <- log1p(-pp) + A - Cval


    if (deriv >= 1) {

      A1 <- ff * (yy - mm) / DD
      dpsi0 <- digamma(yy + ff) - digamma(ff)
      A3 <- ff * (dpsi0 + log(ff / DD) + (mm - yy) / DD)

      A11 <- -mm * ff * (yy + ff) / DD^2
      A13 <- mm * ff * (yy - mm) / DD^2
      dpsi1 <- trigamma(yy + ff) - trigamma(ff)
      A33 <- ff * (dpsi0 + log(ff / DD) + (mm - yy) / DD) +
        ff^2 * dpsi1 + ff * (mm^2 + yy * ff) / DD^2

      l1[ipos, 1L] <- A1 - C1_g[ipos]
      l1[ipos, 2L] <- -pp
      l1[ipos, 3L] <- A3 - C3_g[ipos]

      l2[ipos, 1L] <- A11 - C11_g[ipos]
      l2[ipos, 2L] <- 0
      l2[ipos, 3L] <- A13 - C13_g[ipos]
      l2[ipos, 4L] <- -pp * oo
      l2[ipos, 5L] <- 0
      l2[ipos, 6L] <- A33 - C33_g[ipos]


      if (deriv > 1) {

        A111 <- mm * ff * (yy + ff) * (mm - ff) / DD^3
        A113 <- -mm * ff * (yy + ff) * (mm - ff) / DD^3 - mm * ff^2 / DD^2
        A133 <- mm * ff * (yy - mm) * (mm - ff) / DD^3
        dpsi2 <- psigamma(yy + ff, deriv = 2) - psigamma(ff, deriv = 2)
        B <- dpsi0 + log(ff / DD) + (mm - yy) / DD
        A333 <- ff * B + 3 * ff^2 * dpsi1 + ff^3 * dpsi2 +
          ff * (2 * mm^3 + 3 * mm * ff * yy + ff^2 * yy) / DD^3

        l3[ipos, 1L] <- A111 - C111_g[ipos]
        l3[ipos, 2L] <- 0
        l3[ipos, 3L] <- A113 - C113_g[ipos]
        l3[ipos, 4L] <- 0
        l3[ipos, 5L] <- 0
        l3[ipos, 6L] <- A133 - C133_g[ipos]
        l3[ipos, 7L] <- -pp * oo * (1 - 2 * pp)
        l3[ipos, 8L] <- 0
        l3[ipos, 9L] <- 0
        l3[ipos, 10L] <- A333 - C333_g[ipos]


        if (deriv > 3) {

          C4 <- (mm - ff)^2 - 2 * mm * ff

          A1111 <- -mm * ff * (yy + ff) * C4 / DD^4
          A1113 <- mm * ff * (yy + ff) * C4 / DD^4 +
            mm * ff^2 * (mm - ff) / DD^3
          A1133 <- -mm * ff * (yy + ff) * C4 / DD^4 -
            mm * ff^2 * (3 * mm - ff) / DD^3
          A1333 <- mm * ff * (yy - mm) * C4 / DD^4
          dpsi3 <- psigamma(yy + ff, deriv = 3) - psigamma(ff, deriv = 3)
          A3333 <- ff * B + 7 * ff^2 * dpsi1 + 6 * ff^3 * dpsi2 + ff^4 * dpsi3 +
            ff * (mm^2 * (3 * mm^2 - 2 * mm * ff + ff^2) +
                    yy * ff * (7 * mm^2 + 2 * mm * ff + ff^2)) / DD^4

          l4[ipos, 1L] <- A1111 - C1111_g[ipos]
          l4[ipos, 2L] <- 0
          l4[ipos, 3L] <- A1113 - C1113_g[ipos]
          l4[ipos, 4L] <- 0
          l4[ipos, 5L] <- 0
          l4[ipos, 6L] <- A1133 - C1133_g[ipos]
          l4[ipos, 7L] <- 0
          l4[ipos, 8L] <- 0
          l4[ipos, 9L] <- 0
          l4[ipos, 10L] <- A1333 - C1333_g[ipos]
          l4[ipos, 11L] <- -pp * oo * (1 - 6 * pp + 6 * pp^2)
          l4[ipos, 12L] <- 0
          l4[ipos, 13L] <- 0
          l4[ipos, 14L] <- 0
          l4[ipos, 15L] <- A3333 - C3333_g[ipos]
        }
      }
    }
  }


  ## ================================================================
  ## ZERO COUNTS: ell_0 = log(p) -- exact, trivial, no q-weighting
  ## ================================================================

  if (any(iz)) {

    pp <- p[iz]
    oo <- om[iz]

    l[iz] <- log(pp)

    if (deriv >= 1) {

      l1[iz, 1L] <- 0
      l1[iz, 2L] <- oo            ## d/d eta2 log(p) = 1-p
      l1[iz, 3L] <- 0

      l2[iz, 1L] <- 0
      l2[iz, 2L] <- 0
      l2[iz, 3L] <- 0
      l2[iz, 4L] <- -pp * oo
      l2[iz, 5L] <- 0
      l2[iz, 6L] <- 0

      if (deriv > 1) {

        l3[iz, ] <- 0
        l3[iz, 7L] <- -pp * oo * (1 - 2 * pp)

        if (deriv > 3) {
          l4[iz, ] <- 0
          l4[iz, 11L] <- -pp * oo * (1 - 6 * pp + 6 * pp^2)
        }
      }
    }
  }


  ## ================================================================
  ## EXPECTED HESSIAN (El2)
  ##
  ## (eta1,eta1), (eta1,eta3): fully closed-form (no intractable term).
  ## (eta1,eta2), (eta2,eta3): exactly zero.
  ## (eta2,eta2): exactly -p(1-p).
  ## (eta3,eta3): per current decision, observed curvature l2[,6]
  ## directly (same intractable-E_NB[psi_1] issue as ziNBll's Abar33;
  ## revisit alongside that if/when a fast approximation is settled).
  ## ================================================================

  if (deriv >= 1) {

    Dc_g <- 1 - w_g   ## P(Y>0 | NB) -- the hurdle truncation mass

    Abar11_g <- -mu * phi / D   ## E_NB[A11(Y)]
    ## Abar13 = 0, so the (eta1,eta3) formula below has no bar-term to add

    El2[, 1L] <- (om / Dc_g) * (Abar11_g - w_g * a11_g) - om * C11_g
    El2[, 2L] <- 0
    El2[, 3L] <- (om / Dc_g) * (0 - w_g * a13_g) - om * C13_g
    El2[, 4L] <- -p * om
    El2[, 5L] <- 0
    El2[, 6L] <- l2[, 6L]
  }


  ## ------------------------------------------------
  ## Return
  ## ------------------------------------------------

  list(
    l   = l,
    l1  = l1,
    l2  = l2,
    El2 = El2,
    l3  = l3,
    l4  = l4
  )
}

#' Hurdle negative binomial location-scale model family
#'
#' @description
#' The \code{hurdNBlss} family implements a hurdle negative binomial model in which one linear
#' predictor controls the mean of the (untruncated) negative binomial, a second controls the
#' probability of a zero, and a third controls the overdispersion. Useable only with
#' \code{\link[mgcv]{gam}}, the three linear predictors are specified via a list of formulae.
#' It should be used with care: simply having a large number of zeroes is not an indication that
#' a hurdle (or zero inflation) component is needed.
#'
#' \code{hurdNBll} evaluates the log-likelihood and its derivatives with respect to the three
#' linear predictors. It is used by \code{hurdNBlss} and is exported only to allow external
#' use; it is not usually called directly.
#'
#' Requires non-negative integer count data with at least one count greater than 1.
#'
#' @details
#' \bold{Specifying the model.} \code{gam} is called with a list containing three formulae. The
#' first has the response on the left hand side and the structure of the linear predictor for the
#' mean \eqn{\mu}{mu} of the negative binomial on the right hand side. The second is one sided and
#' specifies the linear predictor for the probability \eqn{p}{p} of a zero. The third is one sided
#' and specifies the linear predictor for the overdispersion \eqn{\phi}{phi}; use \code{~ 1} for
#' a single (scalar) overdispersion parameter.
#'
#' \bold{Parameterisation.} The negative binomial is in its NB2 form,
#' \deqn{f_{NB}(y;\mu,\phi)=\frac{\Gamma(y+\phi)}{\Gamma(\phi)\,y!}\left(\frac{\phi}{\mu+\phi}\right)^{\phi}\left(\frac{\mu}{\mu+\phi}\right)^{y},\quad y=0,1,2,\ldots}{f_NB(y; mu, phi) = Gamma(y + phi) / (Gamma(phi) y!) * (phi / (mu + phi))^phi * (mu / (mu + phi))^y}
#' with mean \eqn{\mu}{mu} and variance \eqn{\mu+\mu^2/\phi}{mu + mu^2 / phi}, so that \eqn{\phi}{phi}
#' is the size parameter: small \eqn{\phi}{phi} means strong overdispersion and \eqn{\phi\to\infty}{phi -> Inf}
#' recovers the Poisson. The three linear predictors are
#' \deqn{\eta_{1}=\log\mu,\qquad \eta_{2}=\log\frac{p}{1-p},\qquad \eta_{3}=\log\phi,}{eta_1 = log(mu), eta_2 = log(p / (1 - p)), eta_3 = log(phi),}
#' i.e. \eqn{\mu=\exp(\eta_1)}{mu = exp(eta_1)}, \eqn{p=1/(1+\exp(-\eta_2))}{p = 1 / (1 + exp(-eta_2))} and
#' \eqn{\phi=\exp(\eta_3)}{phi = exp(eta_3)}. Only identity links are available, because the
#' parameterisation is directly in terms of \eqn{\log\mu}{log(mu)}, \eqn{\mathrm{logit}\,p}{logit(p)} and
#' \eqn{\log\phi}{log(phi)}.
#'
#' \bold{Likelihood.} Here \eqn{p}{p} is exactly the probability of a zero. Positive counts follow a
#' zero-truncated negative binomial:
#' \deqn{P(Y=0)=p,\qquad P(Y=y)=(1-p)\,\frac{f_{NB}(y;\mu,\phi)}{1-f_{NB}(0;\mu,\phi)},\quad y=1,2,\ldots}{P(Y = 0) = p, P(Y = y) = (1 - p) f_NB(y; mu, phi) / (1 - f_NB(0; mu, phi)), y = 1, 2, ...}
#' where \eqn{f_{NB}(0;\mu,\phi)=(\phi/(\mu+\phi))^{\phi}}{f_NB(0; mu, phi) = (phi / (mu + phi))^phi}. The
#' log-likelihood contribution of observation \eqn{i}{i} is
#' \deqn{\ell_i=\log p_i\ \ (y_i=0),\qquad \ell_i=\log(1-p_i)+\log f_{NB}(y_i;\mu_i,\phi_i)-\log\left(1-f_{NB}(0;\mu_i,\phi_i)\right)\ \ (y_i>0).}{l_i = log(p_i) (y_i = 0); l_i = log(1 - p_i) + log f_NB(y_i; mu_i, phi_i) - log(1 - f_NB(0; mu_i, phi_i)) (y_i > 0).}
#' The mean is \eqn{E(Y)=(1-p)\mu/(1-f_{NB}(0;\mu,\phi))}{E(Y) = (1 - p) mu / (1 - f_NB(0; mu, phi))}.
#' Because the zero part depends on \eqn{\eta_2}{eta_2} only, the likelihood separates: the
#' zero-probability coefficients do not interact with the coefficients of the mean and
#' overdispersion predictors (the Hessian is block diagonal), and the mean and overdispersion are
#' estimated from the positive counts alone. Compare \code{\link{ziNBlss}}, in which \eqn{p}{p}
#' is a mixing probability for structural zeros and the zeros couple all three predictors.
#'
#' \bold{Fitted values and prediction.} The fitted values for this family are a three column matrix of the
#' linear predictors \eqn{(\log\mu,\ \mathrm{logit}\,p,\ \log\phi)}{(log(mu), logit(p), log(phi))}.
#' \code{predict.gam} with \code{type = "link"} returns the same three columns, and
#' \code{type = "response"} returns the expected value \eqn{E(Y)}{E(Y)} (with standard errors by the
#' delta method if requested).
#'
#' \bold{Residuals.} \code{residuals} supports \code{"deviance"}, \code{"pearson"}, \code{"scaled.pearson"},
#' \code{"working"} and \code{"response"}. Deviance residuals use the saturated log-likelihood: zero
#' for \eqn{y=0}{y = 0} and \eqn{y=1}{y = 1}, and the zero-truncated Poisson limit
#' (\eqn{\phi\to\infty}{phi -> Inf}, \eqn{p=0}{p = 0}) for \eqn{y>1}{y > 1}.
#'
#' \bold{Null deviance.} The null deviance assumes that a single probability of a zero, a single negative
#' binomial mean and a single overdispersion parameter are estimated (the probability of a zero in closed
#' form, the other two by optimisation).
#'
#' \bold{Identification.} The overdispersion is informed only by the positive counts, and can be weakly
#' identified, particularly when most positive counts equal 1 or the mean is small. In such cases the
#' likelihood can be nearly flat along a direction in which \eqn{\log\mu}{log(mu)} and
#' \eqn{\log\phi}{log(phi)} decrease together (the zero-truncated negative binomial approaches the
#' logarithmic series distribution as \eqn{\phi\to 0}{phi -> 0}), so extremely small estimates of
#' \eqn{\phi}{phi} should be checked. For data with large areas of covariate space over which the response is
#' zero, low order penalties (e.g. \code{s(x, m = 1)}) can help avoid extreme estimates.
#'
#' @section Bias correction of the dispersion parameter (\code{correct.bias}):
#' In \code{mgcv}'s general family framework the unpenalised coefficients of every linear predictor,
#' including the intercept of the overdispersion predictor, are estimated jointly with the penalised
#' coefficients as a penalised-likelihood mode, and the Laplace-approximate REML criterion is used only to choose
#' the smoothing parameters. Marginal-likelihood fitters such as \pkg{glmmTMB} and \pkg{sdmTMB} instead
#' integrate out the smooth coefficients (and, under REML, the mean-predictor coefficients) but not the
#' dispersion coefficients, which makes the dispersion estimate from \code{mgcv} larger (less
#' overdispersed), by an amount that grows with the effective degrees of freedom of the mean
#' predictor. With \code{correct.bias = TRUE} (the default) a one-step correction of the unpenalised
#' coefficients of the overdispersion predictor towards the marginal-likelihood estimate is applied after
#' fitting, in the \code{postproc} step, using \code{n.newton = 0} (only the unpenalised coefficients are
#' shifted). For this family the zero-probability predictor is not corrected, since its likelihood
#' separates from the rest. The coefficients before correction are stored in the fitted object as
#' \code{joint_coef} and the full output of \code{\link{bias_correct}} as \code{bias.correction};
#' \code{Vp}, \code{edf}, \code{sp}, \code{aic} and \code{gcv.ubre} are those of the uncorrected (joint-mode)
#' fit. If the correction fails, the uncorrected fit is returned and the error is stored in \code{bc}.
#' The correction assumes smoothing parameter selection by \code{method = "REML"} (or \code{"ML"}); set
#' \code{correct.bias = FALSE} to obtain the uncorrected fit. See \code{\link{bias_correct}} for details.
#'
#' @param link a list of three character strings specifying the links. Currently only
#'   \code{"identity"} is available, as the parameterisation is directly in terms of
#'   \eqn{\log\mu}{log(mu)}, \eqn{\mathrm{logit}\,p}{logit(p)} and \eqn{\log\phi}{log(phi)}.
#' @param correct.bias logical. If \code{TRUE} (default), apply the one-step marginal-likelihood
#'   correction of the unpenalised coefficients of the overdispersion linear predictor after fitting. See
#'   the section on bias correction and \code{\link{bias_correct}}.
#' @param y response: a vector of non-negative integer counts.
#' @param eta1,eta2,eta3 numeric vectors, of the same length as \code{y}, holding the three linear predictors
#'   \eqn{\log\mu}{log(mu)}, \eqn{\mathrm{logit}\,p}{logit(p)} and \eqn{\log\phi}{log(phi)}.
#' @param deriv integer between 0 and 4 giving the number of derivatives to compute: 0 returns the
#'   log-likelihood only; 1 additionally the first and second derivatives (and the expected Hessian);
#'   2 additionally the third derivatives; 4 additionally the fourth derivatives. (The value 3 is accepted
#'   and returns the same arrays as 2.)
#'
#' @return
#' For \code{hurdNBlss}, an object inheriting from class \code{general.family}.
#'
#' For \code{hurdNBll}, a list with elements
#' \describe{
#'   \item{\code{l}}{vector of log-likelihood contributions.}
#'   \item{\code{l1}}{\eqn{n\times 3}{n x 3} matrix of first derivatives with respect to
#'     \eqn{\eta_1,\eta_2,\eta_3}{eta_1, eta_2, eta_3}.}
#'   \item{\code{l2}}{\eqn{n\times 6}{n x 6} matrix of second derivatives, columns ordered 11, 12, 13, 22, 23, 33.
#'     Columns 12 and 23 are identically zero in this model.}
#'   \item{\code{l3}}{\eqn{n\times 10}{n x 10} matrix of third derivatives, columns ordered
#'     111, 112, 113, 122, 123, 133, 222, 223, 233, 333.}
#'   \item{\code{l4}}{\eqn{n\times 15}{n x 15} matrix of fourth derivatives, columns ordered
#'     1111, 1112, 1113, 1122, 1123, 1133, 1222, 1223, 1233, 1333, 2222, 2223, 2233, 2333, 3333.}
#'   \item{\code{El2}}{\eqn{n\times 6}{n x 6} matrix holding the expected Hessian in the same layout as
#'     \code{l2}. The (\eqn{\eta_3,\eta_3}{eta_3, eta_3}) element is the observed curvature, as the exact
#'     expectation involves an intractable negative binomial term.}
#' }
#' Elements that are not requested by \code{deriv} are \code{NULL}. All derivatives are with respect to the
#' linear predictors \eqn{\eta_1,\eta_2,\eta_3}{eta_1, eta_2, eta_3}.
#'
#' @section WARNINGS:
#' Zero inflated and hurdle models are often over-used. Having many zeroes in the data does not in itself
#' imply that they are needed; having too many zeroes \emph{given the model mean} may.
#'
#' @references
#' Mullahy, J. (1986) Specification and testing of some modified count data models.
#' \emph{Journal of Econometrics} 33, 341-365.
#'
#' Wood, S.N., N. Pya and B. Saefken (2016) Smoothing parameter and model selection for general smooth
#' models. \emph{Journal of the American Statistical Association} 111, 1548-1575.
#' \doi{10.1080/01621459.2016.1180986}
#'
#' @seealso \code{\link{ziNBlss}}, \code{\link{ziPlss}}, \code{\link{bias_correct}},
#'   \code{\link[mgcv]{ziplss}}, \code{\link[mgcv]{gam}}, \code{\link[mgcv]{general.family}}
#'
#' @examples
#' library(mgcv)
#' set.seed(1)
#' n <- 1000
#' x1 <- runif(n); x2 <- runif(n)
#' mu  <- exp(0.5 + sin(2 * pi * x1))        # mean of the untruncated NB
#' p   <- plogis(-0.5 + 2 * (x2 - 0.5))      # P(Y = 0)
#' phi <- 2
#' ## zero-truncated NB draws by inversion on (f0, 1)
#' f0 <- dnbinom(0, mu = mu, size = phi)
#' yp <- qnbinom(runif(n, f0, 1), mu = mu, size = phi)
#' y  <- ifelse(runif(n) < p, 0, yp)
#'
#' \donttest{
#' b <- gam(list(y ~ s(x1), ~ s(x2), ~ 1), family = hurdNBlss(), method = "REML")
#' b$outer.info                 # check convergence
#' summary(b)
#' plot(b, pages = 1)
#' ## coefficients before correction are in b$joint_coef
#' ## b0 <- gam(list(y ~ s(x1), ~ s(x2), ~ 1), family = hurdNBlss(correct.bias = FALSE),
#' ##           method = "REML")
#' }
#' @rdname hurdNBlss
#' @export
hurdNBlss <- function (link = list("identity", "identity", "identity"), correct.bias = TRUE)
{
  if (length(link) != 3)
    stop("hurdNBlss requires 3 links specified as character strings")
  okLinks <- list("identity", "identity", "identity")
  stats <- list()
  param.names <- c("neg. bin. mean", "binary probability", "overdispersion")
  for (i in 1:3) {
    if (link[[i]] %in% okLinks[[i]])
      stats[[i]] <- make.link(link[[i]])
    else stop(link[[i]], " link not available for hurdNBlss")
    fam <- structure(list(link = link[[i]], canonical = "none",
                          linkfun = stats[[i]]$linkfun, mu.eta = stats[[i]]$mu.eta),
                     class = "family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }
  residuals <- function(object,
                        type = c("deviance", "pearson", "scaled.pearson",
                                 "working", "response")) {
    type <- match.arg(type)

    eta1 <- object$fitted.values[, 1]   # log mu
    eta2 <- object$fitted.values[, 2]   # logit p, p = P(y = 0)
    eta3 <- object$fitted.values[, 3]   # log phi
    y    <- object$y

    mu  <- exp(eta1)
    p   <- plogis(eta2)
    phi <- exp(eta3)

    ## NB zero probability and 1 - f0, computed stably for small mu
    lf0  <- phi * (log(phi) - log(phi + mu))
    f0   <- exp(lf0)
    omf0 <- pmax(-expm1(lf0), .Machine$double.eps)

    ## hurdle mean: E[y] = (1 - p) * mu / (1 - f0)
    Ey  <- (1 - p) * mu / omf0
    rsd <- y - Ey
    if (type == "response") return(rsd)

    if (type %in% c("pearson", "scaled.pearson")) {
      ## E[y^2 | y > 0] = (mu + mu^2 * (1 + 1/phi)) / (1 - f0)
      m2 <- (mu + mu^2 * (1 + 1 / phi)) / omf0
      Vy <- pmax((1 - p) * m2 - Ey^2, .Machine$double.eps)
      ## no free scale parameter, so scaled.pearson == pearson
      return(rsd / sqrt(Vy))
    }

    if (type == "working") {
      ## response residual / dE[y]/d eta1
      ## dE/d eta1 = Ey * (1 - mu * phi * f0 / ((phi + mu) * (1 - f0)))
      dE <- Ey * (1 - mu * phi * f0 / ((phi + mu) * omf0))
      return(rsd / dE)
    }

    ## deviance residuals
    ## Saturated log-lik: y = 0 -> p = 1, l = 0.
    ## y = 1 -> the truncated pmf can reach 1, l = 0.
    ## y > 1 -> p = 0 and the zero-truncated Poisson limit (phi -> Inf), with
    ##   g solving g / (1 - exp(-g)) = y  (the same root as the old glo table)
    ls <- numeric(length(y))
    for (v in unique(y[y > 1])) {
      g <- uniroot(function(g) g / (-expm1(-g)) - v,
                   interval = c(v - 1, v), tol = 1e-12)$root
      ls[y == v] <- v * log(g) - g - lgamma(v + 1) - log(-expm1(-g))
    }

    l <- hurdNBll(y, eta1, eta2, eta3, deriv = 0)$l
    sign(rsd) * sqrt(pmax(0, 2 * (ls - l)))
  }
  predict <- function(family, se = FALSE, eta = NULL, y = NULL,
                      X = NULL, beta = NULL, off = NULL, Vb = NULL) {
    if (is.null(eta)) {
      if (is.null(off)) off <- list(0, 0, 0)
      off[[4]] <- 0
      for (i in 1:3) if (is.null(off[[i]])) off[[i]] <- 0
      lpi <- attr(X, "lpi")

      eta1 <- drop(X[, lpi[[1]], drop = FALSE] %*% beta[lpi[[1]]] + off[[1]])  ## log mu
      eta2 <- drop(X[, lpi[[2]], drop = FALSE] %*% beta[lpi[[2]]] + off[[2]])  ## logit p
      eta3 <- drop(X[, lpi[[3]], drop = FALSE] %*% beta[lpi[[3]]] + off[[3]])  ## log phi
    } else {
      se <- FALSE
      eta1 <- eta[, 1]; eta2 <- eta[, 2]; eta3 <- eta[, 3]
    }

    mu  <- exp(eta1)
    p   <- plogis(eta2)
    phi <- exp(eta3)

    ## NB zero probability, computed stably for small mu
    lf0  <- phi * (log(phi) - log(phi + mu))
    f0   <- exp(lf0)
    omf0 <- pmax(-expm1(lf0), .Machine$double.eps)

    ## hurdle mean: E[y] = (1 - p) * mu / (1 - f0)
    Ey <- (1 - p) * mu / omf0
    fv <- list(Ey)

    if (!se) return(fv)

    ## delta method. Derivatives of E[y] w.r.t. the three linear predictors:
    ##   d/d eta1 = Ey * (1 - f0 * phi * mu / ((phi + mu) * (1 - f0)))
    ##   d/d eta2 = -Ey * p
    ##   d/d eta3 = Ey * f0 * A / (1 - f0),
    ##              A = phi * (log(phi / (phi + mu)) + mu / (phi + mu))   (A < 0)
    A  <- phi * (log(phi / (phi + mu)) + mu / (phi + mu))
    d1 <- Ey * (1 - f0 * phi * mu / ((phi + mu) * omf0))
    d2 <- -Ey * p
    d3 <- Ey * f0 * A / omf0

    ## gradient w.r.t. the full coefficient vector, then se = sqrt(diag(G Vb G'))
    G <- matrix(0, nrow(X), ncol(X))
    G[, lpi[[1]]] <- X[, lpi[[1]], drop = FALSE] * d1
    G[, lpi[[2]]] <- X[, lpi[[2]], drop = FALSE] * d2
    G[, lpi[[3]]] <- X[, lpi[[3]], drop = FALSE] * d3
    se.fit <- sqrt(pmax(0, rowSums((G %*% Vb) * G)))

    fv[[2]] <- se.fit
    names(fv) <- c("fit", "se.fit")
    fv
  }
  rd <- function(mu, wt, scale) {
    ## mu: n x 3 matrix of linear predictors (log mean, logit p, log size),
    ## p = P(y = 0) exactly
    if (NCOL(mu) < 3) stop("rd for hurdNBlss needs the 3 linear predictors")
    n   <- nrow(mu)
    m   <- exp(mu[, 1])
    p   <- plogis(mu[, 2])
    phi <- exp(mu[, 3])

    y  <- numeric(n)
    ip <- runif(n) >= p                     # TRUE -> positive count
    np <- sum(ip)
    if (np > 0) {
      m1 <- m[ip]; phi1 <- phi[ip]
      f0 <- dnbinom(0, mu = m1, size = phi1)   # NB P(Y = 0)
      ## inverse CDF restricted to (f0, 1): gives the zero-truncated NB
      u <- runif(np, f0, 1)
      one.eps <- 1 - .Machine$double.eps^0.75
      u[u > one.eps] <- one.eps
      yp <- qnbinom(u, mu = m1, size = phi1)
      y[ip] <- pmax(yp, 1)                  # guard against u == f0 rounding
    }
    y
  }
  postproc <- expression({

    # TODO: it is possible that we want the bias correction after the subsequent deviance calculation but I haven't given this much thought
    # Bias correction for unpenalised coefficients in the higher-order predictors
    if (object$family$correct.bias) {
      bc <- try(bias_correct(object, G, n.newton = 0L))
      if(!inherits(bc, "try-error")) {
        object$joint_coef <- bc$coefficients.joint
        names(object$joint_coef) <- G$term.names
        object$coefficients <- bc$coefficients
        object$bias.correction <- bc
        object$linear.predictors[] <- bc$eta
        object$fitted.values[]     <- bc$eta
      } else {
        object$bc <- bc
      }
    }

    ## 1. Saturated log-likelihood (same as in `residuals`; it is repeated here
    ##    because postproc is evaluated in mgcv's frames, not the constructor's)
    ##    y = 0: p = 1, l = 0.  y = 1: truncated pmf -> 1, l = 0.
    ##    y > 1: zero-truncated Poisson limit, g solving g / (1 - exp(-g)) = y
    ls_hurd <- function(y) {
      ls <- numeric(length(y))
      for (v in unique(y[y > 1])) {
        g <- uniroot(function(g) g / (-expm1(-g)) - v,
                     interval = c(v - 1, v), tol = 1e-12)$root
        ls[y == v] <- v * log(g) - g - lgamma(v + 1) - log(-expm1(-g))
      }
      ls
    }

    ## 2. Intercept-only null model.
    ##    In a hurdle model the likelihood splits into a Bernoulli part for
    ##    y == 0 and a zero-truncated NB part for y > 0, so the null MLE of p
    ##    is closed-form. Only (log mu, log phi) need optimising, which avoids
    ##    the flat p/phi ridge of the ZINB null.
    y_obs <- object$y
    p_obs <- min(max(mean(y_obs == 0), 1e-8), 1 - 1e-8)
    eta_p0 <- qlogis(p_obs)

    hurd_null_obj <- function(params, y, eta_p0) {
      n_obs <- length(y)
      s <- sum(hurdNBll(y, rep(params[1], n_obs), rep(eta_p0, n_obs),
                        rep(params[2], n_obs), deriv = 0)$l)
      if (!is.finite(s)) return(1e10)
      -s
    }

    ypos  <- y_obs[y_obs > 0]
    m_y   <- if (length(ypos)) mean(ypos) else 1
    v_y   <- if (length(ypos) > 1) var(ypos) else 0
    start_mu  <- log(max(0.1, m_y))
    start_phi <- if (is.finite(v_y) && v_y > m_y) log(m_y^2 / (v_y - m_y)) else log(10)

    null_fit <- optim(par = c(start_mu, start_phi), fn = hurd_null_obj,
                      y = y_obs, eta_p0 = eta_p0, method = "Nelder-Mead")
    lnull <- -null_fit$value

    ## 3. Deviances, using the fitted unpenalised log-likelihood
    lsat <- sum(ls_hurd(y_obs))
    lhat <- sum(hurdNBll(y_obs,
                         object$fitted.values[, 1],
                         object$fitted.values[, 2],
                         object$fitted.values[, 3],
                         deriv = 0)$l)

    object$null.deviance <- 2 * (lsat - lnull)
    object$deviance      <- 2 * (lsat - lhat)
  })
  ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL,
                  Hi = NULL, R = NULL, offset = NULL, dH = NULL, db = NULL,
                  deriv = FALSE, nt = 1) {
    gamlss.ncv(X, y, wt, nei, beta, family, llf, H = H, Hi = Hi,
               R = R, offset = offset, dH = dH, db = db, deriv = deriv,
               nt = nt)
  }
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL, D = NULL,
                 eta = NULL, ncv = FALSE, sandwich = FALSE) {
    if (is.null(offset))
      offset <- list(0, 0, 0)
    else offset[[4]] <- 0
    for (i in 1:3) if (is.null(offset[[i]]))
      offset[[i]] <- 0
    jj <- attr(X, "lpi")
    if (is.null(eta)) {
      eta <- X[, jj[[1]], drop = FALSE] %*% coef[jj[[1]]] +
        offset[[1]]
      eta1 <- X[, jj[[2]], drop = FALSE] %*% coef[jj[[2]]] +
        offset[[2]]
      eta2 <- X[, jj[[3]], drop = FALSE] %*% coef[jj[[3]]] +
        offset[[3]]
    }
    else {
      eta2 <- eta[, 3]
      eta1 <- eta[, 2]
      eta <- eta[, 1]
    }
    mu <- family$linfo[[1]]$linkinv(eta)
    p <- family$linfo[[2]]$linkinv(eta1)
    phi <- family$linfo[[3]]$linkinv(eta2)
    zl <- hurdNBll(y, mu, p, phi, deriv)
    # print(paste("dim(zinbll output) = ", Sys.time()))
    # print(lapply(zl, dim))
    if (deriv) {
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta), family$linfo[[2]]$mu.eta(eta1), family$linfo[[3]]$mu.eta(eta2))
      g2 <- cbind(family$linfo[[1]]$d2link(mu), family$linfo[[2]]$d2link(p), family$linfo[[3]]$d2link(phi))
    }
    g3 <- g4 <- matrix(0, nrow(X), 3)
    if (deriv > 1) {
      g3 <- cbind(family$linfo[[1]]$d3link(mu), family$linfo[[2]]$d3link(p), family$linfo[[3]]$d3link(phi))
    }
    if (deriv > 3) {
      g4 <- cbind(family$linfo[[1]]$d4link(mu), family$linfo[[2]]$d4link(p), family$linfo[[3]]$d4link(phi))
    }
    if (deriv) {
      i2 <- family$tri$i2
      i3 <- family$tri$i3
      i4 <- family$tri$i4
      # print(paste("manually set: ", Sys.time()))
      # print(lapply(list(ig1 = ig1, g2 = g2, g3 = g3, g4 = g4, i2 = i2, i3 = i3, i4 = i4), dim))
      de <- gamlss.etamu(zl$l1, zl$l2, zl$l3, zl$l4, ig1,
                         g2, g3, g4, i2, i3, i4, deriv - 1)
      # print(paste("dim(gamlss.etamu output) = ", Sys.time()))
      # print(lapply(de, dim))
      # print(paste("manually set: ", Sys.time()))
      # print(paste0("dim(X) = ", dim(X)))
      # print(paste0("length(jj) = ", length(jj)))
      # print(paste("Columns in X for mu:", length(jj[[1]])))
      # print(paste("Columns in X for p:",  length(jj[[2]])))
      # print(paste("Columns in X for phi:", length(jj[[3]])))
      # print(paste("Total sum of jj tracks:", length(jj[[1]]) + length(jj[[2]]) + length(jj[[3]])))
      ret <- gamlss.gH(X, jj, de$l1, de$l2, i2, l3 = de$l3,
                       i3 = i3, l4 = de$l4, i4 = i4, d1b = d1b, d2b = d2b,
                       deriv = deriv - 1, fh = fh, D = D, sandwich = sandwich)
      if (ncv) {
        ret$l1 <- de$l1
        ret$l2 = de$l2
        ret$l3 = de$l3
      }
    }
    else ret <- list()
    ret$l <- sum(zl$l)
    ret
  }
  sandwich <- function(y, X, coef, wt, family, offset = NULL) {
    ll(y, X, coef, wt, family, offset = NULL, deriv = 1,
       sandwich = TRUE)$lbb
  }
  initialize <- expression({
    n <- rep(1, nobs)
    if (any(y < 0) || !isTRUE(all.equal(y, round(y)))) {
      stop("hurdNBlss requires non-negative integer responses")
    }
    if (max(y) < 2) stop("hurdNBlss needs some counts above 1 to identify the truncated NB")
    if (is.null(start)) {
      use.unscaled <- !is.null(attr(E, "use.unscaled"))
      jj <- attr(x, "lpi")

      start <- (function(y, x, E, jj, use.unscaled) {
        ## ---- helpers ----------------------------------------------------
        ## (penalised) least squares, optionally weighted, same two branches as before
        pen_fit <- function(xx, ee, yy, ww = NULL) {
          if (!is.null(ww)) {
            ww <- pmax(ww, 1e-8); ww <- ww / mean(ww)
            sw <- sqrt(ww); xx <- xx * sw; yy <- yy * sw
          }
          if (use.unscaled) {
            cf <- qr.coef(qr(rbind(xx, ee)), c(yy, rep(0, nrow(ee))))
            cf[!is.finite(cf)] <- 0
            cf
          } else pen.reg(xx, ee, yy)
        }
        lf0 <- function(mu, phi) -phi * log1p(mu / phi)            # log NB P(Y = 0)
        ## zero-truncated NB log-likelihood (mean mu, size phi)
        tnb_l <- function(yy, mu, phi)
          dnbinom(yy, size = phi, mu = mu, log = TRUE) -
          log(pmax(-expm1(lf0(mu, phi)), 1e-300))
        ## solve mu / (1 - f0(mu, phi)) = m for mu by vectorised bisection.
        ## The truncated mean increases in mu, tends to 1 as mu -> 0, and is >= mu.
        inv_trunc <- function(m, phi) {
          m   <- pmax(m, 1 + 1e-6)
          phi <- rep_len(phi, length(m))
          lo  <- rep(log(1e-8), length(m)); hi <- log(m)
          for (k in 1:40) {
            mid <- (lo + hi) / 2; mu <- exp(mid)
            h   <- mu / pmax(-expm1(lf0(mu, phi)), 1e-300)
            up  <- h > m
            hi[up] <- mid[up]; lo[!up] <- mid[!up]
          }
          exp((lo + hi) / 2)
        }

        st  <- rep(0, ncol(x))
        pos <- y > 0
        yp  <- y[pos]; npos <- length(yp); mbar <- mean(yp)

        ## ---- 1. global truncated-NB MLE of (log mu, log phi), positives only ----
        par0 <- c(log(inv_trunc(mbar, 1)), 0)
        gfit <- tryCatch(
          optim(par0, function(pr) {
            v <- -sum(tnb_l(yp, exp(pr[1]), exp(pr[2])))
            if (is.finite(v)) v else 1e10
          }, method = "L-BFGS-B",
          lower = c(-12, -5), upper = c(log(10 * max(yp) + 10), 12))$par,
          error = function(e) par0)
        if (any(!is.finite(gfit))) gfit <- par0
        phi_g <- exp(gfit[2])

        ## ---- 2. mean predictor: quasi-likelihood IRLS for the truncated mean m(x) ----
        x1p <- x[pos, jj[[1]], drop = FALSE]
        e1  <- E[, jj[[1]], drop = FALSE]
        eta_m <- rep(log(max(mbar, 1 + 1e-6)), npos)
        if (npos >= 10) {
          for (k in 1:3) {
            m_i <- exp(eta_m)
            zz  <- eta_m + (yp - m_i) / m_i
            cf  <- pen_fit(x1p, e1, zz, m_i / (1 + m_i / phi_g))
            eta_m <- pmin(pmax(drop(x1p %*% cf), 1e-6), log(max(yp)) + 1)
          }
        }

        ## ---- 3. refine phi given the fitted means, then convert m -> mu ----
        mu_i <- inv_trunc(exp(eta_m), phi_g)
        lp <- tryCatch(
          optimize(function(l) -sum(tnb_l(yp, mu_i, exp(l))), c(-5, 12))$minimum,
          error = function(e) gfit[2])
        mu_i <- inv_trunc(exp(eta_m), exp(lp))
        st[jj[[1]]] <- pen_fit(x1p, e1, log(mu_i))

        ## ---- 4. dispersion predictor: constant start at the refined phi ----
        x3 <- x[, jj[[3]], drop = FALSE]
        e3 <- E[, jj[[3]], drop = FALSE]
        st[jj[[3]]] <- pen_fit(x3, e3, rep(lp, nrow(x3)))

        ## ---- 5. zero probability: binomial IRLS on all rows, p = P(y = 0) ----
        x2  <- x[, jj[[2]], drop = FALSE]
        e2  <- E[, jj[[2]], drop = FALSE]
        z0  <- as.numeric(y == 0)
        eta_p <- rep(qlogis(min(max(mean(z0), 1e-4), 1 - 1e-4)), nrow(x))
        for (k in 1:3) {
          p_i <- plogis(eta_p); w_i <- pmax(p_i * (1 - p_i), 1e-6)
          cf  <- pen_fit(x2, e2, eta_p + (z0 - p_i) / w_i, w_i)
          eta_p <- pmin(pmax(drop(x2 %*% cf), -10), 10)
        }
        st[jj[[2]]] <- cf
        st
      })(y, x, E, jj, use.unscaled)
    }
  })
  structure(list(family = "hurdNBlss", ll = ll, link = paste(link), correct.bias = correct.bias,
                 nlp = 3, ncv = ncv, sandwich = sandwich, tri = trind.generator(3),
                 initialize = initialize, postproc = postproc, residuals = residuals,
                 rd = rd, predict = predict, linfo = stats, d2link = 1,
                 d3link = 1, d4link = 1, ls = 1, available.derivs = 2),
            class = c("general.family", "extended.family", "family"))
}

## Mixture Version #############################################################

# log-likelihood function and derivatives

#' @rdname ziNBlss
#' @export
ziNBll <- function(y, eta1, eta2, eta3, deriv = 0) {

  ## ================================================================
  ## Zero-inflated negative binomial log likelihood
  ##
  ## Parameterisation:
  ##
  ##   mu  = exp(eta1)
  ##   p   = plogis(eta2)
  ##   phi = exp(eta3)
  ##
  ## NB2:
  ##
  ##   Var(Y | NB) = mu + mu^2 / phi
  ##
  ## The function returns derivatives with respect to
  ## eta1, eta2 and eta3.
  ##
  ## deriv:
  ##   0 : likelihood only
  ##   1 : likelihood, first and second derivatives
  ##   2 : additionally third derivatives
  ##   3/4 : additionally fourth derivatives
  ##
  ## No Taylor-jet or automatic-differentiation calculations are used.
  ##
  ## ----------------------------------------------------------------
  ## NOTE on El2 (expected Hessian):
  ##
  ## Unlike l1..l4 (observed derivatives, which necessarily depend on
  ## the realised y), the expected Hessian El2 is a property of the
  ## model at (mu, phi, p) only and must NOT depend on which y was
  ## actually observed at a given row. It is therefore computed once,
  ## for every row, from the shared "starred" (y=0) NB quantities
  ## (Astar, A1s, A3s, A11s, A13s) together with q = p/P0, rather than
  ## being branched by iz/ipos.
  ##
  ## The (eta3,eta3) entry still uses the observed curvature l2[,6]
  ## as a placeholder (mirroring mgcv::negbin's treatment of the NB
  ## dispersion curvature), because the exact expected value requires
  ## the intractable E_NB[psi_1(Y+phi)] term. Swap that placeholder
  ## out once the approximation for that term is finalised.
  ## ================================================================


  ## ------------------------------------------------
  ## Input checks
  ## ------------------------------------------------

  n <- length(y)

  if (length(eta1) != n ||
      length(eta2) != n ||
      length(eta3) != n) {
    stop("y, eta1, eta2 and eta3 must have equal lengths")
  }

  if (length(deriv) != 1L ||
      !is.finite(deriv) ||
      deriv != as.integer(deriv) ||
      deriv < 0 || deriv > 4) {
    stop("deriv must be an integer between 0 and 4")
  }

  if (any(!is.finite(y)) ||
      any(y < 0) ||
      any(y != floor(y))) {
    stop("y must contain finite non-negative integers")
  }

  deriv <- as.integer(deriv)


  ## ------------------------------------------------
  ## Stable inverse links
  ## ------------------------------------------------

  ## mu and phi require exponentiation. The clipping is only to
  ## prevent floating-point overflow/underflow in the parameter
  ## transformations themselves.

  lo <- log(.Machine$double.xmin)
  hi <- log(.Machine$double.xmax)

  e1 <- pmin(pmax(eta1, lo), hi)
  e3 <- pmin(pmax(eta3, lo), hi)

  mu  <- exp(e1)
  phi <- exp(e3)

  ## Stable logistic transformation
  p <- numeric(n)

  ip <- eta2 >= 0
  p[ip]  <- 1 / (1 + exp(-eta2[ip]))

  ep <- exp(eta2[!ip])
  p[!ip] <- ep / (1 + ep)

  om <- 1 - p


  ## ------------------------------------------------
  ## Allocate outputs
  ## ------------------------------------------------

  l <- numeric(n)

  l1 <- l2 <- l3 <- l4 <- El2 <- NULL

  if (deriv >= 1) {
    l1 <- matrix(0, n, 3L)
    l2 <- matrix(0, n, 6L)
    El2 <- matrix(0, n, 6L)
  }

  if (deriv > 1) {
    l3 <- matrix(0, n, 10L)
  }

  if (deriv > 3) {
    l4 <- matrix(0, n, 15L)
  }


  ## ------------------------------------------------
  ## Common NB quantities
  ## ------------------------------------------------

  D <- mu + phi

  log_phi_over_D <- log(phi) - log(D)


  ## ------------------------------------------------
  ## Row-wise "starred" (y = 0) quantities, needed for El2.
  ##
  ## These depend only on mu, phi, p -- NOT on the observed y --
  ## and are therefore computed for every row, up front, once.
  ## ------------------------------------------------

  if (deriv >= 1) {

    ## A* = A(0; eta1, eta3)
    Astar_g <- phi * log_phi_over_D

    A1s_g  <- -mu * phi / D
    A3s_g  <- phi * (log_phi_over_D + mu / D)
    A11s_g <- -mu * phi^2 / D^2
    A13s_g <- -mu^2 * phi / D^2

    ## P0 = P(Y = 0) under the full ZINB model
    f0_g <- exp(Astar_g)
    P0_g <- p + om * f0_g

    ## q = P(structural zero | Y = 0) = p / P0 = logit^{-1}(eta2 - A*)
    xg <- eta2 - Astar_g

    q_g <- numeric(n)
    iqg <- xg >= 0
    q_g[iqg]  <- 1 / (1 + exp(-xg[iqg]))
    exg <- exp(xg[!iqg])
    q_g[!iqg] <- exg / (1 + exg)

    qm_g <- 1 - q_g
    q1_g <- q_g * qm_g
  }


  ## ================================================================
  ## POSITIVE COUNTS
  ## ================================================================

  iz <- y == 0
  ipos <- !iz

  if (any(ipos)) {

    yy  <- y[ipos]
    mm  <- mu[ipos]
    pp  <- p[ipos]
    oo  <- om[ipos]
    ff  <- phi[ipos]
    DD  <- D[ipos]

    ## --------------------------------------------------------------
    ## NB log likelihood A(y; mu, phi)
    ##
    ## A =
    ## log Gamma(y+phi) - log Gamma(phi)
    ## + phi log(phi/(mu+phi))
    ## + y log(mu/(mu+phi))
    ## --------------------------------------------------------------

    A <- lgamma(yy + ff) -
      lgamma(ff) -
      lgamma(yy + 1) +
      ff * log(ff / DD) +
      yy * log(mm / DD)

    l[ipos] <- log1p(-pp) + A


    if (deriv >= 1) {

      ## ------------------------------------------------------------
      ## First derivatives of A
      ## ------------------------------------------------------------

      A1 <- ff * (yy - mm) / DD

      dpsi0 <- digamma(yy + ff) - digamma(ff)

      A3 <- ff * (
        dpsi0 +
          log(ff / DD) +
          (mm - yy) / DD
      )


      ## ------------------------------------------------------------
      ## Second derivatives
      ## ------------------------------------------------------------

      A11 <- -mm * ff * (yy + ff) / DD^2

      A13 <- mm * ff * (yy - mm) / DD^2

      dpsi1 <- trigamma(yy + ff) - trigamma(ff)

      A33 <- ff * (
        dpsi0 +
          log(ff / DD) +
          (mm - yy) / DD
      ) +
        ff^2 * dpsi1 +
        ff * (mm^2 + yy * ff) / DD^2


      ## ------------------------------------------------------------
      ## First derivatives of log likelihood
      ## ------------------------------------------------------------

      l1[ipos, 1L] <- A1

      ## d/d eta2 log(1-p) = -p
      l1[ipos, 2L] <- -pp

      l1[ipos, 3L] <- A3


      ## ------------------------------------------------------------
      ## Second derivatives
      ##
      ## ordering:
      ## 11, 12, 13, 22, 23, 33
      ## ------------------------------------------------------------

      l2[ipos, 1L] <- A11
      l2[ipos, 2L] <- 0
      l2[ipos, 3L] <- A13

      l2[ipos, 4L] <- -pp * oo
      l2[ipos, 5L] <- 0

      l2[ipos, 6L] <- A33


      if (deriv > 1) {

        ## ----------------------------------------------------------
        ## Third derivatives of A
        ## ----------------------------------------------------------

        A111 <-
          mm * ff * (yy + ff) * (mm - ff) / DD^3

        A113 <-
          -mm * ff * (yy + ff) * (mm - ff) / DD^3 -
          mm * ff^2 / DD^2

        A133 <-
          mm * ff * (yy - mm) * (mm - ff) / DD^3

        dpsi2 <- psigamma(yy + ff, deriv = 2) -
          psigamma(ff, deriv = 2)

        B <- dpsi0 +
          log(ff / DD) +
          (mm - yy) / DD

        A333 <-
          ff * B +
          3 * ff^2 * dpsi1 +
          ff^3 * dpsi2 +
          ff * (
            2 * mm^3 +
              3 * mm * ff * yy +
              ff^2 * yy
          ) / DD^3


        ## ----------------------------------------------------------
        ## Third derivatives of positive-count likelihood
        ##
        ## ordering:
        ## 111,112,113,122,123,133,222,223,233,333
        ## ----------------------------------------------------------

        l3[ipos, 1L] <- A111
        l3[ipos, 2L] <- 0
        l3[ipos, 3L] <- A113

        l3[ipos, 4L] <- 0
        l3[ipos, 5L] <- 0
        l3[ipos, 6L] <- A133

        l3[ipos, 7L] <- -pp * oo * (1 - 2 * pp)
        l3[ipos, 8L] <- 0
        l3[ipos, 9L] <- 0

        l3[ipos, 10L] <- A333


        if (deriv > 3) {

          ## --------------------------------------------------------
          ## Fourth derivatives of A
          ## --------------------------------------------------------

          C4 <- (mm - ff)^2 - 2 * mm * ff

          A1111 <-
            -mm * ff * (yy + ff) * C4 / DD^4

          A1113 <-
            mm * ff * (yy + ff) * C4 / DD^4 +
            mm * ff^2 * (mm - ff) / DD^3

          A1133 <-
            -mm * ff * (yy + ff) * C4 / DD^4 -
            mm * ff^2 * (3 * mm - ff) / DD^3

          A1333 <-
            mm * ff * (yy - mm) * C4 / DD^4

          dpsi3 <- psigamma(yy + ff, deriv = 3) -
            psigamma(ff, deriv = 3)

          A3333 <-
            ff * B +
            7 * ff^2 * dpsi1 +
            6 * ff^3 * dpsi2 +
            ff^4 * dpsi3 +
            ff * (
              mm^2 * (
                3 * mm^2 -
                  2 * mm * ff +
                  ff^2
              ) +
                yy * ff * (
                  7 * mm^2 +
                    2 * mm * ff +
                    ff^2
                )
            ) / DD^4


          ## --------------------------------------------------------
          ## Fourth derivatives of positive-count likelihood
          ##
          ## ordering:
          ## 1111,1112,1113,1122,1123,1133,
          ## 1222,1223,1233,1333,
          ## 2222,2223,2233,2333,3333
          ## --------------------------------------------------------

          l4[ipos, 1L] <- A1111
          l4[ipos, 2L] <- 0
          l4[ipos, 3L] <- A1113

          l4[ipos, 4L] <- 0
          l4[ipos, 5L] <- 0
          l4[ipos, 6L] <- A1133

          l4[ipos, 7L] <- 0
          l4[ipos, 8L] <- 0
          l4[ipos, 9L] <- 0
          l4[ipos, 10L] <- A1333

          l4[ipos, 11L] <-
            -pp * oo * (1 - 6 * pp + 6 * pp^2)

          l4[ipos, 12L] <- 0
          l4[ipos, 13L] <- 0
          l4[ipos, 14L] <- 0
          l4[ipos, 15L] <- A3333
        }
      }
    }
  }


  ## ================================================================
  ## ZERO COUNTS
  ## ================================================================

  if (any(iz)) {

    mm <- mu[iz]
    pp <- p[iz]
    oo <- om[iz]
    ff <- phi[iz]
    DD <- D[iz]

    ## --------------------------------------------------------------
    ## NB probability of zero
    ##
    ## A* = phi log(phi/(mu+phi))
    ## --------------------------------------------------------------

    Ast <- ff * (log(ff) - log(DD))

    ## --------------------------------------------------------------
    ## Stable zero-inflated log likelihood
    ##
    ## l0 = A* - softplus(eta2)
    ##      + softplus(eta2 - A*)
    ##
    ## --------------------------------------------------------------

    x <- eta2[iz] - Ast

    ## Stable softplus
    softplus_x <- pmax(x, 0) + log1p(exp(-abs(x)))
    softplus_eta2 <-
      pmax(eta2[iz], 0) +
      log1p(exp(-abs(eta2[iz])))

    l[iz] <-
      Ast -
      softplus_eta2 +
      softplus_x


    if (deriv >= 1) {

      ## ------------------------------------------------------------
      ## Starred NB derivatives
      ## ------------------------------------------------------------

      A1 <- -mm * ff / DD

      A3 <- ff * (
        log(ff / DD) +
          mm / DD
      )

      A11 <- -mm * ff^2 / DD^2

      A13 <- -mm^2 * ff / DD^2

      ## A33 = the general A33 formula evaluated at y=0
      ## (dpsi0 = dpsi1 = 0 there, so the psi terms drop out).

      A33 <- ff * (
        log(ff / DD) +
          mm / DD
      ) +
        ff * mm^2 / DD^2

      ## ------------------------------------------------------------
      ## Posterior structural-zero probability
      ##
      ## q = p / [p + (1-p) exp(A*)]
      ##   = logit^{-1}(eta2 - A*)
      ##
      ## Already computed for every row above as q_g; subset it here
      ## for use in the observed zero-count derivatives.
      ## ------------------------------------------------------------

      q  <- q_g[iz]
      qm <- qm_g[iz]
      q1 <- q1_g[iz]
      q2 <- q1 * (1 - 2 * q)
      q3 <- q1 * (1 - 6 * q + 6 * q^2)

      ## ------------------------------------------------------------
      ## First derivatives
      ## ------------------------------------------------------------

      l1[iz, 1L] <- qm * A1
      l1[iz, 2L] <- q - pp
      l1[iz, 3L] <- qm * A3


      ## ------------------------------------------------------------
      ## Second derivatives
      ##
      ## 11,12,13,22,23,33
      ## ------------------------------------------------------------

      l2[iz, 1L] <- qm * A11 + q1 * A1^2

      l2[iz, 2L] <- -q1 * A1

      l2[iz, 3L] <- qm * A13 + q1 * A1 * A3

      l2[iz, 4L] <- q1 - pp * oo

      l2[iz, 5L] <- -q1 * A3

      l2[iz, 6L] <- qm * A33 + q1 * A3^2


      if (deriv > 1) {

        ## ----------------------------------------------------------
        ## Additional starred NB third derivatives
        ## ----------------------------------------------------------

        A111 <-
          mm * ff^2 * (mm - ff) / DD^3

        A113 <-
          -mm * ff^2 * (mm - ff) / DD^3 -
          mm * ff^2 / DD^2

        A133 <-
          -mm^2 * ff * (mm - ff) / DD^3

        A333 <-
          ff * (
            log(ff / DD) +
              mm / DD
          ) +
          2 * ff * mm^3 / DD^3


        ## ----------------------------------------------------------
        ## Third derivatives of zero-count likelihood
        ##
        ## 111,112,113,122,123,133,
        ## 222,223,233,333
        ## ----------------------------------------------------------

        l3[iz, 1L] <-
          qm * A111 +
          3 * q1 * A1 * A11 -
          q2 * A1^3

        l3[iz, 2L] <-
          q2 * A1^2 -
          q1 * A11

        l3[iz, 3L] <-
          qm * A113 +
          q1 * (
            A3 * A11 +
              2 * A1 * A13
          ) -
          q2 * A1^2 * A3

        l3[iz, 4L] <-
          -q2 * A1

        l3[iz, 5L] <-
          q2 * A1 * A3 -
          q1 * A13

        l3[iz, 6L] <-
          qm * A133 +
          q1 * (
            A1 * A33 +
              2 * A3 * A13
          ) -
          q2 * A1 * A3^2

        l3[iz, 7L] <-
          q2 -
          pp * oo * (1 - 2 * pp)

        l3[iz, 8L] <-
          -q2 * A3

        l3[iz, 9L] <-
          -q1 * A33 +
          q2 * A3^2

        l3[iz, 10L] <-
          qm * A333 +
          3 * q1 * A3 * A33 -
          q2 * A3^3


        if (deriv > 3) {

          ## --------------------------------------------------------
          ## Additional starred NB fourth derivatives
          ## --------------------------------------------------------

          C4 <- (mm - ff)^2 - 2 * mm * ff

          A1111 <-
            -mm * ff^2 * C4 / DD^4

          A1113 <-
            mm * ff^2 * C4 / DD^4 +
            mm * ff^2 * (mm - ff) / DD^3

          A1133 <-
            -mm * ff^2 * C4 / DD^4 -
            mm * ff^2 * (3 * mm - ff) / DD^3

          A1333 <-
            -mm^2 * ff * C4 / DD^4

          A3333 <-
            ff * (
              log(ff / DD) +
                mm / DD
            ) +
            ff * mm^2 * (
              3 * mm^2 -
                2 * mm * ff +
                ff^2
            ) / DD^4


          ## --------------------------------------------------------
          ## Fourth derivatives of zero-count likelihood
          ##
          ## 1111,1112,1113,1122,1123,1133,
          ## 1222,1223,1233,1333,
          ## 2222,2223,2233,2333,3333
          ## --------------------------------------------------------

          l4[iz, 1L] <-
            qm * A1111 +
            q1 * (
              4 * A1 * A111 +
                3 * A11^2
            ) -
            6 * q2 * A1^2 * A11 +
            q3 * A1^4

          l4[iz, 2L] <-
            -q1 * A111 +
            3 * q2 * A1 * A11 -
            q3 * A1^3

          l4[iz, 3L] <-
            qm * A1113 +
            q1 * (
              3 * A1 * A113 +
                A3 * A111 +
                3 * A11 * A13
            ) -
            3 * q2 * (
              A1^2 * A13 +
                A1 * A3 * A11
            ) +
            q3 * A1^3 * A3

          l4[iz, 4L] <-
            q3 * A1^2 -
            q2 * A11

          l4[iz, 5L] <-
            -q1 * A113 +
            q2 * (
              2 * A1 * A13 +
                A3 * A11
            ) -
            q3 * A1^2 * A3

          l4[iz, 6L] <-
            qm * A1133 +
            q1 * (
              2 * A1 * A133 +
                2 * A3 * A113 +
                A11 * A33 +
                2 * A13^2
            ) -
            q2 * (
              A1^2 * A33 +
                4 * A1 * A3 * A13 +
                A3^2 * A11
            ) +
            q3 * A1^2 * A3^2

          l4[iz, 7L] <-
            -q3 * A1

          l4[iz, 8L] <-
            -q2 * A13 +
            q3 * A1 * A3

          l4[iz, 9L] <-
            -q1 * A133 +
            q2 * (
              A1 * A33 +
                2 * A3 * A13
            ) -
            q3 * A1 * A3^2

          l4[iz, 10L] <-
            qm * A1333 +
            q1 * (
              A1 * A333 +
                3 * A3 * A133 +
                3 * A13 * A33
            ) -
            3 * q2 * (
              A1 * A3 * A33 +
                A13 * A3^2
            ) +
            q3 * A1 * A3^3

          l4[iz, 11L] <-
            q3 -
            pp * oo * (1 - 6 * pp + 6 * pp^2)

          l4[iz, 12L] <-
            -q3 * A3

          l4[iz, 13L] <-
            q3 * A3^2 -
            q2 * A33

          l4[iz, 14L] <-
            -q1 * A333 +
            3 * q2 * A3 * A33 -
            q3 * A3^3

          l4[iz, 15L] <-
            qm * A3333 +
            q1 * (
              4 * A3 * A333 +
                3 * A33^2
            ) -
            6 * q2 * A3^2 * A33 +
            q3 * A3^4
        }
      }
    }
  }


  ## ================================================================
  ## EXPECTED HESSIAN (El2)
  ##
  ## Assembled once, uniformly across ALL rows, using only the
  ## row-wise quantities mu, phi, p, D, Astar_g, A1s_g, A3s_g, q_g
  ## computed above -- never branched by iz/ipos, since the expected
  ## Hessian cannot depend on which y was actually observed.
  ##
  ## Closed forms (derived from D*(1-q) = (1-p)*exp(Astar), which
  ## cancels the exp(Astar)*A_ij* correction terms exactly):
  ##
  ##   H_11 = (1-p)*Abar11 + P0*q*(1-q)*A1s^2      Abar11 = -mu*phi/(mu+phi)
  ##   H_12 = -P0*q*(1-q)*A1s
  ##   H_13 =  P0*q*(1-q)*A1s*A3s                   (Abar13 = 0 drops out)
  ##   H_22 =  P0*q*(1-q) - p*(1-p)
  ##   H_23 = -P0*q*(1-q)*A3s
  ##   H_33 = (1-p)*Abar33 + P0*q*(1-q)*A3s^2       (placeholder below)
  ##
  ## where P0 = P(Y=0) under the full ZINB model (matches D from the
  ## write-up, but named P0_g here to avoid clashing with mu+phi).
  ## ================================================================

  if (deriv >= 1) {

    El2[, 1L] <- om * (-mu * phi / D) + P0_g * q1_g * A1s_g^2

    El2[, 2L] <- -P0_g * q1_g * A1s_g

    El2[, 3L] <- P0_g * q1_g * A1s_g * A3s_g

    El2[, 4L] <- P0_g * q1_g - p * om

    El2[, 5L] <- -P0_g * q1_g * A3s_g

    ## --------------------------------------------------------------
    ## (eta3, eta3) placeholder: observed curvature, following
    ## mgcv::negbin's treatment of the NB dispersion term. Abar33
    ## requires the intractable E_NB[psi_1(Y+phi)] term -- once an
    ## approximation for that is finalised, replace this line with
    ## the full closed form:
    ##
    ##   (1 - p) * Abar33 + P0_g * q1_g * A3s_g^2
    ## --------------------------------------------------------------

    El2[, 6L] <- pmin(l2[,6], -1e-6)
  }


  ## ------------------------------------------------
  ## Return
  ## ------------------------------------------------

  list(
    l   = l,
    l1  = l1,
    l2  = l2,
    El2 = El2,
    l3  = l3,
    l4  = l4
  )
}

# family constructor

#' Zero-inflated (mixture) negative binomial location-scale model family
#'
#' @description
#' The \code{ziNBlss} family implements a zero-inflated negative binomial model in which one linear
#' predictor controls the mean of the negative binomial, a second controls the probability that an
#' observation is a structural zero, and a third controls the overdispersion. Useable only with
#' \code{\link[mgcv]{gam}}, the three linear predictors are specified via a list of formulae. It should be used
#' with care: simply having a large number of zeroes is not an indication of zero inflation.
#'
#' \code{ziNBll} evaluates the log-likelihood and its derivatives with respect to the three linear
#' predictors. It is used by \code{ziNBlss} and is exported only to allow external use; it is not usually
#' called directly.
#'
#' Requires non-negative integer count data.
#'
#' @details
#' \bold{Specifying the model.} \code{gam} is called with a list containing three formulae. The first has the
#' response on the left hand side and the structure of the linear predictor for the negative binomial mean
#' \eqn{\mu}{mu} on the right hand side. The second is one sided and specifies the linear predictor for the
#' probability \eqn{p}{p} of a structural zero. The third is one sided and specifies the linear predictor for the
#' overdispersion \eqn{\phi}{phi}; use \code{~ 1} for a single (scalar) overdispersion parameter.
#'
#' \bold{Parameterisation.} The negative binomial is in its NB2 form,
#' \deqn{f_{NB}(y;\mu,\phi)=\frac{\Gamma(y+\phi)}{\Gamma(\phi)\,y!}\left(\frac{\phi}{\mu+\phi}\right)^{\phi}\left(\frac{\mu}{\mu+\phi}\right)^{y},\quad y=0,1,2,\ldots}{f_NB(y; mu, phi) = Gamma(y + phi) / (Gamma(phi) y!) * (phi / (mu + phi))^phi * (mu / (mu + phi))^y}
#' with mean \eqn{\mu}{mu} and variance \eqn{\mu+\mu^2/\phi}{mu + mu^2 / phi}; \eqn{\phi}{phi} is the size
#' parameter, small values meaning strong overdispersion and \eqn{\phi\to\infty}{phi -> Inf} recovering the
#' Poisson. The three linear predictors are
#' \deqn{\eta_{1}=\log\mu,\qquad \eta_{2}=\log\frac{p}{1-p},\qquad \eta_{3}=\log\phi,}{eta_1 = log(mu), eta_2 = log(p / (1 - p)), eta_3 = log(phi),}
#' i.e. \eqn{\mu=\exp(\eta_1)}{mu = exp(eta_1)}, \eqn{p=1/(1+\exp(-\eta_2))}{p = 1 / (1 + exp(-eta_2))} and
#' \eqn{\phi=\exp(\eta_3)}{phi = exp(eta_3)}. Only identity links are available, because the parameterisation
#' is directly in terms of \eqn{\log\mu}{log(mu)}, \eqn{\mathrm{logit}\,p}{logit(p)} and
#' \eqn{\log\phi}{log(phi)}.
#'
#' \bold{Likelihood.} Here \eqn{p}{p} is a mixing probability: with probability \eqn{p}{p} an observation is a
#' structural zero, and otherwise it is drawn from the negative binomial, so that zeros can arise either way:
#' \deqn{P(Y=0)=p+(1-p)\,f_{NB}(0;\mu,\phi),\qquad P(Y=y)=(1-p)\,f_{NB}(y;\mu,\phi),\quad y=1,2,\ldots}{P(Y = 0) = p + (1 - p) f_NB(0; mu, phi), P(Y = y) = (1 - p) f_NB(y; mu, phi), y = 1, 2, ...}
#' where \eqn{f_{NB}(0;\mu,\phi)=(\phi/(\mu+\phi))^{\phi}}{f_NB(0; mu, phi) = (phi / (mu + phi))^phi}. The
#' log-likelihood contribution of observation \eqn{i}{i} is
#' \deqn{\ell_i=\log\left(p_i+(1-p_i)f_{NB}(0;\mu_i,\phi_i)\right)\ \ (y_i=0),\qquad \ell_i=\log(1-p_i)+\log f_{NB}(y_i;\mu_i,\phi_i)\ \ (y_i>0).}{l_i = log(p_i + (1 - p_i) f_NB(0; mu_i, phi_i)) (y_i = 0); l_i = log(1 - p_i) + log f_NB(y_i; mu_i, phi_i) (y_i > 0).}
#' The mean and variance are \eqn{E(Y)=(1-p)\mu}{E(Y) = (1 - p) mu} and
#' \eqn{\mathrm{Var}(Y)=(1-p)\mu\,(1+\mu/\phi+p\mu)}{Var(Y) = (1 - p) mu (1 + mu / phi + p mu)}. In the
#' derivatives for zero counts the posterior probability that a zero is structural,
#' \eqn{q=p/(p+(1-p)f_{NB}(0;\mu,\phi))}{q = p / (p + (1 - p) f_NB(0; mu, phi))}, appears throughout. Unlike the
#' hurdle model of \code{\link{hurdNBlss}}, zero counts couple all three linear predictors, so the Hessian
#' is not block diagonal. For the same reason \eqn{p}{p} and \eqn{\phi}{phi} both govern the probability of a zero,
#' and with little information they can be only weakly identified separately.
#'
#' \bold{Fitted values and prediction.} The fitted values for this family are a three column matrix of the linear
#' predictors \eqn{(\log\mu,\ \mathrm{logit}\,p,\ \log\phi)}{(log(mu), logit(p), log(phi))}.
#' \code{predict.gam} with \code{type = "link"} returns the same three columns, and
#' \code{type = "response"} returns \eqn{E(Y)}{E(Y)} (with standard errors by the delta method if requested; these
#' depend on the mean and zero-inflation predictors only).
#'
#' \bold{Residuals.} \code{residuals} supports \code{"deviance"}, \code{"pearson"} and \code{"response"}. Deviance
#' residuals use the saturated log-likelihood: zero for \eqn{y=0}{y = 0}, and for \eqn{y>0}{y > 0} the Poisson
#' limit (\eqn{p=0}{p = 0}, \eqn{\mu=y}{mu = y}, \eqn{\phi\to\infty}{phi -> Inf}),
#' \eqn{y\log y-y-\log y!}{y log(y) - y - log(y!)}.
#'
#' \bold{Null deviance.} The null deviance assumes that a single mean, a single probability of a structural zero and a
#' single overdispersion parameter are estimated.
#'
#' @section Bias correction of the zero-inflation and dispersion parameters (\code{correct.bias}):
#' In \code{mgcv}'s general family framework the unpenalised coefficients of every linear predictor (for example
#' the intercepts, and the linear terms of smooths) are estimated jointly with the penalised coefficients as a
#' penalised-likelihood mode, and the Laplace-approximate REML criterion is used only to choose the smoothing
#' parameters. Marginal-likelihood fitters such as \pkg{glmmTMB} and \pkg{sdmTMB} integrate out the smooth
#' coefficients and (under REML) the mean-predictor coefficients, but not the unpenalised coefficients of the
#' zero-inflation and dispersion predictors. For the dispersion this makes the \code{mgcv} estimate larger (less
#' overdispersed), by an amount that grows with the effective degrees of freedom of the mean predictor. With
#' \code{correct.bias = TRUE} (the default) a one-step correction of the unpenalised coefficients of the
#' zero-inflation and dispersion predictors towards the marginal-likelihood estimate is applied after fitting, in
#' the \code{postproc} step, using \code{n.newton = 0} (only the unpenalised coefficients are shifted, jointly,
#' since they are correlated). The coefficients before correction are stored in the fitted object as
#' \code{joint_coef} and the full output of \code{\link{bias_correct}} as \code{bias.correction}; \code{Vp},
#' \code{edf}, \code{sp}, \code{aic} and \code{gcv.ubre} are those of the uncorrected (joint-mode) fit. If the
#' correction fails, the uncorrected fit is returned and the error is stored in \code{bc}. The correction assumes
#' smoothing parameter selection by \code{method = "REML"} (or \code{"ML"}); set \code{correct.bias = FALSE} to
#' obtain the uncorrected fit. See \code{\link{bias_correct}} for details and limitations.
#'
#' @param link a list of three character strings specifying the links. Currently only
#'   \code{"identity"} is available, as the parameterisation is directly in terms of
#'   \eqn{\log\mu}{log(mu)}, \eqn{\mathrm{logit}\,p}{logit(p)} and \eqn{\log\phi}{log(phi)}.
#' @param correct.bias logical. If \code{TRUE} (default), apply the one-step marginal-likelihood correction of the
#'   unpenalised coefficients of the zero-inflation and overdispersion linear predictors after fitting. See the
#'   section on bias correction and \code{\link{bias_correct}}.
#' @param y response: a vector of non-negative integer counts.
#' @param eta1,eta2,eta3 numeric vectors, of the same length as \code{y}, holding the three linear predictors
#'   \eqn{\log\mu}{log(mu)}, \eqn{\mathrm{logit}\,p}{logit(p)} and \eqn{\log\phi}{log(phi)}.
#' @param deriv integer between 0 and 4 giving the number of derivatives to compute: 0 returns the
#'   log-likelihood only; 1 additionally the first and second derivatives (and the expected Hessian);
#'   2 additionally the third derivatives; 4 additionally the fourth derivatives. (The value 3 is accepted
#'   and returns the same arrays as 2.)
#'
#' @return
#' For \code{ziNBlss}, an object inheriting from class \code{general.family}.
#'
#' For \code{ziNBll}, a list with elements
#' \describe{
#'   \item{\code{l}}{vector of log-likelihood contributions.}
#'   \item{\code{l1}}{\eqn{n\times 3}{n x 3} matrix of first derivatives with respect to
#'     \eqn{\eta_1,\eta_2,\eta_3}{eta_1, eta_2, eta_3}.}
#'   \item{\code{l2}}{\eqn{n\times 6}{n x 6} matrix of second derivatives, columns ordered 11, 12, 13, 22, 23, 33.}
#'   \item{\code{l3}}{\eqn{n\times 10}{n x 10} matrix of third derivatives, columns ordered
#'     111, 112, 113, 122, 123, 133, 222, 223, 233, 333.}
#'   \item{\code{l4}}{\eqn{n\times 15}{n x 15} matrix of fourth derivatives, columns ordered
#'     1111, 1112, 1113, 1122, 1123, 1133, 1222, 1223, 1233, 1333, 2222, 2223, 2233, 2333, 3333.}
#'   \item{\code{El2}}{\eqn{n\times 6}{n x 6} matrix holding the expected Hessian in the same layout as
#'     \code{l2}. The (\eqn{\eta_3,\eta_3}{eta_3, eta_3}) element is the observed curvature (bounded away from
#'     zero), as the exact expectation involves an intractable negative binomial term.}
#' }
#' Elements that are not requested by \code{deriv} are \code{NULL}. All derivatives are with respect to the linear
#' predictors \eqn{\eta_1,\eta_2,\eta_3}{eta_1, eta_2, eta_3}.
#'
#' @section WARNINGS:
#' Zero inflated models are often over-used. Having many zeroes in the data does not in itself imply zero
#' inflation. Having too many zeroes \emph{given the model mean} may imply zero inflation.
#'
#' @references
#' Lambert, D. (1992) Zero-inflated Poisson regression, with an application to defects in manufacturing.
#' \emph{Technometrics} 34, 1-14.
#'
#' Wood, S.N., N. Pya and B. Saefken (2016) Smoothing parameter and model selection for general smooth
#' models. \emph{Journal of the American Statistical Association} 111, 1548-1575.
#' \doi{10.1080/01621459.2016.1180986}
#'
#' @seealso \code{\link{hurdNBlss}}, \code{\link{ziPlss}}, \code{\link{bias_correct}},
#'   \code{\link[mgcv]{ziplss}}, \code{\link[mgcv]{gam}}, \code{\link[mgcv]{general.family}}
#'
#' @examples
#' library(mgcv)
#' set.seed(2)
#' n <- 1000
#' x1 <- runif(n); x2 <- runif(n)
#' mu  <- exp(0.5 + sin(2 * pi * x1))
#' p   <- plogis(-1 + 2 * (x2 - 0.5))       # P(structural zero)
#' phi <- 2
#' y <- rnbinom(n, mu = mu, size = phi)
#' y[runif(n) < p] <- 0
#'
#' \donttest{
#' b <- gam(list(y ~ s(x1), ~ s(x2), ~ 1), family = ziNBlss(), method = "REML")
#' b$outer.info
#' summary(b)
#' plot(b, pages = 1)
#' ## uncorrected (joint-mode) fit for comparison
#' b0 <- gam(list(y ~ s(x1), ~ s(x2), ~ 1), family = ziNBlss(correct.bias = FALSE),
#'           method = "REML")
#' cbind(joint = coef(b0), corrected = coef(b))[c("(Intercept).1", "(Intercept).2"), ]
#' }
#' @rdname ziNBlss
#' @export
ziNBlss <- function (link = list("identity", "identity", "identity"), correct.bias = TRUE)
{
  if (length(link) != 3)
    stop("ziNBlss requires 3 links specified as character strings")
  okLinks <- list("identity", "identity", "identity")
  stats <- list()
  param.names <- c("neg. bin. mean", "binary probability", "overdispersion")
  for (i in 1:3) {
    if (link[[i]] %in% okLinks[[i]])
      stats[[i]] <- make.link(link[[i]])
    else stop(link[[i]], " link not available for ziNBlss")
    fam <- structure(list(link = link[[i]], canonical = "none",
                          linkfun = stats[[i]]$linkfun, mu.eta = stats[[i]]$mu.eta),
                     class = "family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }
  residuals <- function(object, type = c("deviance", "pearson", "scaled.pearson",
                                         "working", "response")) {
    type <- match.arg(type)

    eta1 <- object$fitted.values[, 1]   # log mu
    eta2 <- object$fitted.values[, 2]   # logit p (structural zero prob)
    eta3 <- object$fitted.values[, 3]   # log phi (NB size)
    y    <- object$y

    mu  <- exp(eta1)
    p   <- plogis(eta2)
    phi <- exp(eta3)

    Ey  <- (1 - p) * mu
    rsd <- y - Ey

    if (type == "response") return(rsd)

    if (type == "pearson") {
      Vy <- (1 - p) * mu * (1 + mu / phi + p * mu)
      return(rsd / sqrt(Vy))
    }

    ## deviance residuals
    ## Saturated log-lik: y == 0 -> p = 1 gives l = 0.
    ## y > 0 -> p = 0, mu = y, phi -> Inf (the Poisson limit), the supremum
    ## of the NB pmf at y, giving y*log(y) - y - lgamma(y + 1).
    ls <- ifelse(y > 0, y * log(pmax(y, 1)) - y - lgamma(y + 1), 0)

    l   <- ziNBll(y, eta1, eta2, eta3, deriv = 0)$l
    dev <- pmax(0, 2 * (ls - l))
    sign(rsd) * sqrt(dev)
  }
  predict <- function(family, se = FALSE, eta = NULL, y = NULL,
                      X = NULL, beta = NULL, off = NULL, Vb = NULL) {
    if (is.null(eta)) {
      if (is.null(off)) off <- list(0, 0, 0)
      off[[4]] <- 0
      for (i in 1:3) if (is.null(off[[i]])) off[[i]] <- 0
      lpi <- attr(X, "lpi")
      X1 <- X[, lpi[[1]], drop = FALSE]   ## mu block
      X2 <- X[, lpi[[2]], drop = FALSE]   ## p block

      eta1 <- drop(X1 %*% beta[lpi[[1]]] + off[[1]])   ## linear predictor -> mu
      eta2 <- drop(X2 %*% beta[lpi[[2]]] + off[[2]])   ## linear predictor -> p

      mu <- exp(eta1)
      p  <- plogis(eta2)

      if (se) {
        v1  <- drop(pmax(0, rowSums((X1 %*% Vb[lpi[[1]], lpi[[1]]]) * X1)))  ## Var(eta1)
        v2  <- drop(pmax(0, rowSums((X2 %*% Vb[lpi[[2]], lpi[[2]]]) * X2)))  ## Var(eta2)
        v12 <- drop(rowSums((X1 %*% Vb[lpi[[1]], lpi[[2]]]) * X2))           ## Cov(eta1,eta2)
      }
    } else {
      se <- FALSE
      eta1 <- eta[, 1]
      eta2 <- eta[, 2]
      mu <- exp(eta1)
      p  <- plogis(eta2)
    }

    fv <- list((1 - p) * mu)

    if (!se) {
      return(fv)
    } else {
      ## delta method: d(fv)/d(eta1) = (1-p)*mu ; d(fv)/d(eta2) = -mu*p*(1-p)
      dfv.deta1 <- (1 - p) * mu
      dfv.deta2 <- -mu * p * (1 - p)

      se.fit <- sqrt(pmax(0,
                          dfv.deta1^2 * v1 + dfv.deta2^2 * v2 + 2 * dfv.deta1 * dfv.deta2 * v12
      ))

      fv[[2]] <- se.fit
      names(fv) <- c("fit", "se.fit")
      return(fv)
    }
  }
  rd <- function(mu, wt, scale) {
    ## mu: n x 3 matrix of linear predictors (log mean, logit p, log size),
    ## the same convention the ziplss rd used (mu[,1] was log lambda)
    if (NCOL(mu) < 3) stop("rd for ziNBlss needs the 3 linear predictors")
    n   <- nrow(mu)
    m   <- exp(mu[, 1])
    p   <- plogis(mu[, 2])
    phi <- exp(mu[, 3])
    y   <- rnbinom(n, mu = m, size = phi)   # NB draw
    y[runif(n) < p] <- 0                    # structural zeros
    y
  }
  postproc <- expression({

    # TODO: it is possible that we want the bias correction after the subsequent deviance calculation but I haven't given this much thought
    # Bias correction for unpenalised coefficients in the higher-order predictors
    if (object$family$correct.bias) {
      bc <- try(bias_correct(object, G, n.newton = 0L))
      if(!inherits(bc, "try-error")) {
        object$joint_coef <- bc$coefficients.joint
        names(object$joint_coef) <- G$term.names
        object$coefficients <- bc$coefficients
        object$bias.correction <- bc
        object$linear.predictors[] <- bc$eta
        object$fitted.values[]     <- bc$eta
      } else {
        object$bc <- bc
      }
    }

    # 1. Saturated Log-Likelihood for Negative Binomial Space
    # A saturated discrete model fits perfectly when mu_i = y_i, p_i = 0, and phi_i -> infinity.
    # At that mathematical limit, the log-likelihood reduces exactly to 0 for all counts.
    ls <- function(y) {
      rep(0, length(y))
    }

    # 2. Intercept-Only Null Model Profile Function for ZINB
    # This evaluates the log-likelihood of a completely flat, 3-intercept null baseline model
    zinb_null_obj <- function(params, y) {
      # params contains: [eta1_mu, eta2_p, eta3_phi]
      n_obs <- length(y)

      # Map to response parameters via links (note the updates here might increase comp. time a bit)
      mu_null  <- rep(params[1], n_obs)
      p_null   <- rep(params[2], n_obs) # FLAG LINK FN
      phi_null <- rep(params[3], n_obs)

      # Call your standard zinbll function (deriv = 0)
      res <- ziNBll(y, mu_null, p_null, phi_null, deriv = 0)

      # Return negative log-likelihood for minimization
      return(-sum(res$l))
    }

    # 3. Find the Intercept-Only Null Model Likelihood
    # Derive safe, realistic starting values from global data attributes
    m_y <- mean(object$y)
    v_y <- var(object$y)
    p_obs <- sum(object$y == 0) / length(object$y)

    start_mu  <- log(max(0.1, m_y))
    start_p   <- log(max(0.01, min(0.99, p_obs)) / (1 - max(0.01, min(0.99, p_obs)))) # FLAG LINK FN
    start_phi <- if (v_y > m_y) log((m_y^2) / (v_y - m_y)) else log(1)

    # Execute a minor Nelder-Mead optimization to find the accurate null profile match
    null_fit <- optim(par = c(start_mu, start_p, start_phi),
                      fn = zinb_null_obj,
                      y = object$y,
                      method = "Nelder-Mead")

    lnull <- -null_fit$value

    # 4. Populate Deviance Metrics into the final Gam Object
    # Saturated baseline sum is 0, so deviance matches standard log-likelihood gaps
    ls <- function(y) ifelse(y > 0, y * log(pmax(y, 1)) - y - lgamma(y + 1), 0)

    object$null.deviance <- 2 * (sum(ls(object$y)) - lnull)
    # object$deviance      <- 2 * (sum(ls(object$y)) - object$lw)
    lhat <- sum(ziNBll(object$y,
                       object$fitted.values[, 1],
                       object$fitted.values[, 2],
                       object$fitted.values[, 3],
                       deriv = 0)$l)
    object$deviance <- 2 * (sum(ls(object$y)) - lhat)

  })
  ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL,
                  Hi = NULL, R = NULL, offset = NULL, dH = NULL, db = NULL,
                  deriv = FALSE, nt = 1) {
    gamlss.ncv(X, y, wt, nei, beta, family, llf, H = H, Hi = Hi,
               R = R, offset = offset, dH = dH, db = db, deriv = deriv,
               nt = nt)
  }
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL, D = NULL,
                 eta = NULL, ncv = FALSE, sandwich = FALSE) {
    if (is.null(offset))
      offset <- list(0, 0, 0)
    else offset[[4]] <- 0
    for (i in 1:3) if (is.null(offset[[i]]))
      offset[[i]] <- 0
    jj <- attr(X, "lpi")
    if (is.null(eta)) {
      eta <- X[, jj[[1]], drop = FALSE] %*% coef[jj[[1]]] +
        offset[[1]]
      eta1 <- X[, jj[[2]], drop = FALSE] %*% coef[jj[[2]]] +
        offset[[2]]
      eta2 <- X[, jj[[3]], drop = FALSE] %*% coef[jj[[3]]] +
        offset[[3]]
    }
    else {
      eta2 <- eta[, 3]
      eta1 <- eta[, 2]
      eta <- eta[, 1]
    }
    mu <- family$linfo[[1]]$linkinv(eta)
    p <- family$linfo[[2]]$linkinv(eta1)
    phi <- family$linfo[[3]]$linkinv(eta2)
    zl <- ziNBll(y, mu, p, phi, deriv)
    # print(paste("dim(zinbll output) = ", Sys.time()))
    # print(lapply(zl, dim))
    if (deriv) {
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta), family$linfo[[2]]$mu.eta(eta1), family$linfo[[3]]$mu.eta(eta2))
      g2 <- cbind(family$linfo[[1]]$d2link(mu), family$linfo[[2]]$d2link(p), family$linfo[[3]]$d2link(phi))
    }
    g3 <- g4 <- matrix(0, nrow(X), 3)
    if (deriv > 1) {
      g3 <- cbind(family$linfo[[1]]$d3link(mu), family$linfo[[2]]$d3link(p), family$linfo[[3]]$d3link(phi))
    }
    if (deriv > 3) {
      g4 <- cbind(family$linfo[[1]]$d4link(mu), family$linfo[[2]]$d4link(p), family$linfo[[3]]$d4link(phi))
    }
    if (deriv) {
      i2 <- family$tri$i2
      i3 <- family$tri$i3
      i4 <- family$tri$i4
      # print(paste("manually set: ", Sys.time()))
      # print(lapply(list(ig1 = ig1, g2 = g2, g3 = g3, g4 = g4, i2 = i2, i3 = i3, i4 = i4), dim))
      de <- gamlss.etamu(zl$l1, zl$l2, zl$l3, zl$l4, ig1,
                         g2, g3, g4, i2, i3, i4, deriv - 1)
      # print(paste("dim(gamlss.etamu output) = ", Sys.time()))
      # print(lapply(de, dim))
      # print(paste("manually set: ", Sys.time()))
      # print(paste0("dim(X) = ", dim(X)))
      # print(paste0("length(jj) = ", length(jj)))
      # print(paste("Columns in X for mu:", length(jj[[1]])))
      # print(paste("Columns in X for p:",  length(jj[[2]])))
      # print(paste("Columns in X for phi:", length(jj[[3]])))
      # print(paste("Total sum of jj tracks:", length(jj[[1]]) + length(jj[[2]]) + length(jj[[3]])))
      ret <- gamlss.gH(X, jj, de$l1, de$l2, i2, l3 = de$l3,
                       i3 = i3, l4 = de$l4, i4 = i4, d1b = d1b, d2b = d2b,
                       deriv = deriv - 1, fh = fh, D = D, sandwich = sandwich)
      if (ncv) {
        ret$l1 <- de$l1
        ret$l2 = de$l2
        ret$l3 = de$l3
      }
    }
    else ret <- list()
    ret$l <- sum(zl$l)
    ret
  }
  sandwich <- function(y, X, coef, wt, family, offset = NULL) {
    ll(y, X, coef, wt, family, offset = NULL, deriv = 1,
       sandwich = TRUE)$lbb
  }
  initialize <- expression({
    n <- rep(1, nobs)
    if (any(y < 0) || !isTRUE(all.equal(y, round(y)))) {
      stop("Non-negative integer response required with ziNBlss")
    }
    if (min(y) == 0 && max(y) == 1) stop("Using ziNBlss for binary data makes no sense")
    if (is.null(start)) {
      use.unscaled <- !is.null(attr(E, "use.unscaled"))
      jj <- attr(x, "lpi")

      start <- (function(y, x, E, jj, use.unscaled, n_em = 3L) {
        ## ---- helpers (local, so nothing leaks into mgcv's frame) ----
        pen_fit <- function(xx, ee, yy, ww = NULL) {
          if (!is.null(ww)) {
            ww <- pmax(ww, 1e-8); ww <- ww / mean(ww)
            sw <- sqrt(ww); xx <- xx * sw; yy <- yy * sw
          }
          if (use.unscaled) {
            cf <- qr.coef(qr(rbind(xx, ee)), c(yy, rep(0, nrow(ee))))
            cf[!is.finite(cf)] <- 0
            cf
          } else pen.reg(xx, ee, yy)
        }
        ## factorial-moment estimate of phi (same estimator as before)
        mom_phi <- function(mu, p) {
          num <- sum((1 - p) * mu^2)
          den <- sum(y * (y - 1) - (1 - p) * mu^2)
          if (is.finite(num) && is.finite(den) && num > 0 && den > 0)
            min(max(num / den, 1e-2), 1e6) else 1
        }
        lo_mu <- log(1e-3); hi_mu <- log(max(y) + 1)

        x1 <- x[, jj[[1]], drop = FALSE]; e1 <- E[, jj[[1]], drop = FALSE]
        x2 <- x[, jj[[2]], drop = FALSE]; e2 <- E[, jj[[2]], drop = FALSE]
        x3 <- x[, jj[[3]], drop = FALSE]; e3 <- E[, jj[[3]], drop = FALSE]
        nn  <- length(y)
        pos <- y > 0

        ## ---- initial values ----
        ## mean: log y on the positive rows only (zeros carry no usable log-scale information)
        cf1 <- if (sum(pos) >= 5) pen_fit(x1[pos, , drop = FALSE], e1, log(y[pos]))
        else pen_fit(x1, e1, rep(log(mean(y[pos])), nn))
        eta1 <- pmin(pmax(drop(x1 %*% cf1), lo_mu), hi_mu)
        mu   <- exp(eta1)

        ## p: constant start from the excess of zeros over the NB zeros at phi = 1
        f0bar <- mean(1 / (1 + mu))
        pbar  <- min(max((mean(!pos) - f0bar) / (1 - f0bar), 0.01), 0.99)
        cf2   <- pen_fit(x2, e2, rep(qlogis(pbar), nn))
        eta2  <- pmin(pmax(drop(x2 %*% cf2), -10), 10)

        phi <- mom_phi(mu, plogis(eta2))

        ## ---- a few EM-style sweeps (no likelihood evaluations) ----
        for (k in seq_len(n_em)) {
          p  <- plogis(eta2)
          f0 <- exp(-phi * log1p(mu / phi))                 # NB P(Y = 0)
          z  <- ifelse(pos, 0, p / (p + (1 - p) * f0))      # P(structural | y)

          ## p: IRLS step on soft labels z
          w2   <- pmax(p * (1 - p), 1e-6)
          cf2  <- pen_fit(x2, e2, eta2 + (z - p) / w2, w2)
          eta2 <- pmin(pmax(drop(x2 %*% cf2), -10), 10)

          ## mu: NB IRLS step, rows weighted by P(NB-generated)
          w1   <- (1 - z) * mu / (1 + mu / phi)
          cf1  <- pen_fit(x1, e1, eta1 + (y - mu) / mu, w1)
          eta1 <- pmin(pmax(drop(x1 %*% cf1), lo_mu), hi_mu)
          mu   <- exp(eta1)

          phi <- mom_phi(mu, plogis(eta2))
        }

        st <- rep(0, ncol(x))
        st[jj[[1]]] <- cf1
        st[jj[[2]]] <- cf2
        st[jj[[3]]] <- pen_fit(x3, e3, rep(log(phi), nrow(x3)))
        st
      })(y, x, E, jj, use.unscaled)
    }
  })
  structure(list(family = "ziNBlss", ll = ll, link = paste(link), correct.bias = correct.bias,
                 nlp = 3, ncv = ncv, sandwich = sandwich, tri = trind.generator(3),
                 initialize = initialize, postproc = postproc, residuals = residuals,
                 rd = rd, predict = predict, linfo = stats, d2link = 1,
                 d3link = 1, d4link = 1, ls = 1, available.derivs = 2),
            class = c("general.family", "extended.family", "family"))
}

################################################################################
#### Zero-inflated Poisson Model (Mixture) #####################################
################################################################################

# log-likelihood function and derivatives

#' @rdname ziPlss
#' @export
ziPll <- function(y, eta1, eta2, deriv = 0) {
  l1 <- El2 <- l2 <- l3 <- l4 <- NULL
  n <- length(y)
  zind <- y == 0

  # Evaluate positive counts scale using canonical link helpers
  slam <- ldlam(eta1, y, deriv = deriv)
  sp   <- ldp(eta2, deriv = deriv)

  lam <- slam$lam
  p <- sp$p

  # 1. Total Sample Log-Likelihood Vector
  l <- rep(0, n)
  l[zind]  <- log(p[zind] + (1 - p[zind]) * exp(-lam[zind]))
  l[!zind] <- log(1 - p[!zind]) + y[!zind] * eta1[!zind] - lam[!zind] - lgamma(y[!zind] + 1)

  if (deriv) {
    l1 <- matrix(0, n, 2)
    l2 <- matrix(0, n, 3)
    El2 <- matrix(0, n, 3)

    # Inject positive count baselines
    l1[!zind, 1] <- slam$l1[!zind]
    l1[!zind, 2] <- sp$l1[!zind]
    l2[!zind, 1] <- slam$l2[!zind]
    l2[!zind, 2] <- 0
    l2[!zind, 3] <- sp$l2[!zind]

    # --- EVALUATE ZERO OBSERVATIONS (MIXTURE CALC) ---
    A <- exp(l[zind]) # Equal to: p + (1-p)*exp(-lam)
    E <- exp(-lam[zind])
    O <- 1 - p[zind]

    p1 <- p[zind] * O
    p2 <- p1 * (1 - 2 * p[zind])

    E1 <- -lam[zind] * E
    E2 <- (lam[zind]^2 - lam[zind]) * E

    A1 <- O * E1
    A2 <- p1 * (1 - E)
    A11 <- O * E2
    A22 <- p2 * (1 - E)
    A12 <- -p1 * E1

    # 1st order arrays for zeros
    l1[zind, 1] <- A1 / A
    l1[zind, 2] <- A2 / A

    # 2nd order arrays for zeros (Hessian entries)
    l2[zind, 1] <- A11 / A - (l1[zind, 1])^2
    l2[zind, 2] <- A12 / A - (l1[zind, 1] * l1[zind, 2])
    l2[zind, 3] <- A22 / A - (l1[zind, 2])^2

    # --- EXPECTED HESSIAN MATRIX (El2) OVER FULL SAMPLE ---
    # Weighting the zero vs conditional non-zero matrices via expected probability A_full
    # A_full <- exp(l)
    A_full <- p + (1 - p) * exp(-lam)

    # Re-calculate zero matrix elements across all observation space for expected limits
    E_f <- exp(-lam)
    O_f <- 1 - p
    p1_f <- p * O_f
    p2_f <- p1_f * (1 - 2 * p)
    E1_f <- -lam * E_f
    E2_f <- (lam^2 - lam) * E_f

    A1_f <- O_f * E1_f
    A2_f <- p1_f * (1 - E_f)
    A11_f <- O_f * E2_f
    A22_f <- p1_f * (1 - E_f) # fixed bug: p2_f * (1 - E_f)
    A22_f <- p2_f * (1 - E_f)
    A12_f <- -p1_f * E1_f

    l2_z_f1 <- A11_f / A_full - (A1_f / A_full)^2
    l2_z_f2 <- A12_f / A_full - (A1_f * A2_f / A_full^2)
    l2_z_f3 <- A22_f / A_full - (A2_f / A_full)^2

    El2[, 1] <- A_full * l2_z_f1 + (1 - A_full) * (-lam)
    # El2[, 2] <- A_full * l2_z_f2 + (1 - A_full) * 0
    El2[, 2] <- A_full * l2_z_f2
    El2[, 3] <- A_full * l2_z_f3 + (1 - A_full) * (-p * O_f)
  }

  if (deriv > 1) {
    l3 <- matrix(0, n, 4)
    l3[!zind, 1] <- slam$l3[!zind]
    l3[!zind, 4] <- sp$l3[!zind]

    # Higher order chain rule terms for zeros
    p3 <- p1 * (1 - 6 * p[zind] + 6 * p[zind]^2)
    E3 <- (-lam[zind]^3 + 3 * lam[zind]^2 - lam[zind]) * E

    A111 <- O * E3
    A222 <- p3 * (1 - E)
    A112 <- -p1 * E2
    A122 <- -p2 * E1

    l3[zind, 1] <- A111 / A - 3 * A11 * A1 / A^2 + 2 * A1^3 / A^3
    l3[zind, 2] <- A112 / A - (A11 * A2 + 2 * A12 * A1) / A^2 + 2 * A1^2 * A2 / A^3
    l3[zind, 3] <- A122 / A - (A22 * A1 + 2 * A12 * A2) / A^2 + 2 * A1 * A2^2 / A^3
    l3[zind, 4] <- A222 / A - 3 * A22 * A2 / A^2 + 2 * A2^3 / A^3
  }

  if (deriv > 3) {
    l4 <- matrix(0, n, 5)
    l4[!zind, 1] <- slam$l4[!zind]
    l4[!zind, 5] <- sp$l4[!zind]

    p4 <- p1 * (1 - 14 * p[zind] + 36 * p[zind]^2 - 24 * p[zind]^3)
    E4 <- (lam[zind]^4 - 6 * lam[zind]^3 + 7 * lam[zind]^2 - lam[zind]) * E

    A1111 <- O * E4
    A2222 <- p4 * (1 - E)
    A1112 <- -p1 * E3
    A1122 <- -p2 * E2
    A1222 <- -p3 * E1

    l4[zind, 1] <- A1111 / A - (4 * A111 * A1 + 3 * A11^2) / A^2 + 12 * A11 * A1^2 / A^3 - 6 * A1^4 / A^4
    l4[zind, 2] <- A1112 / A - (A111 * A2 + 3 * A112 * A1 + 3 * A11 * A12) / A^2 + 6 * (A11 * A1 * A2 + A12 * A1^2) / A^3 - 6 * A1^3 * A2 / A^4
    l4[zind, 3] <- A1122 / A - (A112 * A2 + A122 * A1 + A11 * A22 + 2 * A12^2) / A^2 + 2 * (A11 * A2^2 + A22 * A1^2 + 4 * A12 * A1 * A2) / A^3 - 6 * A1^2 * A2^2 / A^4
    l4[zind, 4] <- A1222 / A - (A222 * A1 + 3 * A122 * A2 + 3 * A22 * A12) / A^2 + 6 * (A22 * A1 * A2 + A12 * A2^2) / A^3 - 6 * A1 * A2^3 / A^4
    l4[zind, 5] <- A2222 / A - (4 * A222 * A2 + 3 * A22^2) / A^2 + 12 * A22 * A2^2 / A^3 - 6 * A2^4 / A^4
  }

  list(l = l, l1 = l1, l2 = l2, l3 = l3, l4 = l4, El2 = El2)
}

# additional derivative functions for each linear predictor:

#' Derivative helpers for the zero-inflated Poisson likelihood
#'
#' @description
#' Internal helper functions used by \code{\link{ziPll}}. \code{ldlam} gives the log-likelihood derivatives of the
#' Poisson part, \eqn{y\eta_1-e^{\eta_1}}{y eta_1 - exp(eta_1)}, with respect to \eqn{\eta_1=\log\lambda}{eta_1 = log(lambda)}.
#' \code{ldp} gives the derivatives of \eqn{\log(1-p)}{log(1 - p)} with respect to
#' \eqn{\eta_2=\mathrm{logit}\,p}{eta_2 = logit(p)}, which is the part of the log-likelihood that does not depend on
#' the count for positive observations. These are the building blocks of the positive-count terms of the
#' zero-inflated Poisson log-likelihood.
#'
#' @details
#' \code{ldlam} returns \eqn{l_1=y-\lambda}{l1 = y - lambda} and \eqn{l_2=l_3=l_4=-\lambda}{l2 = l3 = l4 = -lambda}, where
#' \eqn{\lambda=\exp(\eta_1)}{lambda = exp(eta_1)}, with \eqn{\eta_1}{eta_1} bounded above to avoid overflow.
#' \code{ldp} returns \eqn{l_1=-p}{l1 = -p}, \eqn{l_2=-p(1-p)}{l2 = -p (1 - p)},
#' \eqn{l_3=-p(1-p)(1-2p)}{l3 = -p (1 - p) (1 - 2p)} and \eqn{l_4=-p(1-p)(1-6p+6p^2)}{l4 = -p (1 - p) (1 - 6p + 6p^2)},
#' where \eqn{p=1/(1+\exp(-\eta_2))}{p = 1 / (1 + exp(-eta_2))}, with \eqn{\eta_2}{eta_2} bounded to avoid overflow.
#'
#' @param eta1 linear predictor for \eqn{\log\lambda}{log(lambda)}.
#' @param eta2 linear predictor for \eqn{\mathrm{logit}\,p}{logit(p)}.
#' @param y response counts.
#' @param deriv number of derivatives to compute: third derivatives are returned if \code{deriv > 1} and fourth
#'   derivatives if \code{deriv > 3}.
#'
#' @return
#' For \code{ldlam}, a list with elements \code{l1}, \code{l2}, \code{l3}, \code{l4} (derivatives with respect to
#' \eqn{\eta_1}{eta_1}; \code{l3} and \code{l4} are \code{NULL} unless requested) and \code{lam}, the Poisson mean.
#'
#' For \code{ldp}, a list with elements \code{l1}, \code{l2}, \code{l3}, \code{l4} (derivatives with respect to
#' \eqn{\eta_2}{eta_2}) and \code{p}, the probability of a structural zero.
#'
#' @seealso \code{\link{ziPlss}}
#' @keywords internal
#' @rdname ziP-helpers
ldlam <- function(eta1, y, deriv = 4) {
  # Bound eta1 to avoid exponential overflow using mgcv limits
  hi <- log(.Machine$double.xmax) / 3
  eta1 <- pmin(eta1, hi)

  lam <- exp(eta1)
  n <- length(eta1)

  l1 <- y - lam
  l2 <- -lam

  l3 <- l4 <- NULL
  if (deriv > 1) {
    l3 <- -lam
  }
  if (deriv > 3) {
    l4 <- -lam
  }

  list(l1 = l1, l2 = l2, l3 = l3, l4 = l4, lam = lam)
}

#' @rdname ziP-helpers
#' @keywords internal
ldp <- function(eta2, deriv = 4) {
  # Stabilise eta2 boundaries to avoid NaN in exp(eta2)
  hi <- log(.Machine$double.xmax) / 3
  lo <- log(.Machine$double.eps) / 3
  eta2 <- pmax(pmin(eta2, hi), lo)

  p <- exp(eta2) / (1 + exp(eta2)) # FLAG LINK FN
  n <- length(eta2)

  l1 <- -p
  l2 <- -p * (1 - p)

  l3 <- l4 <- NULL
  if (deriv > 1) {
    l3 <- -p * (1 - p) * (1 - 2 * p)
  }
  if (deriv > 3) {
    l4 <- -p * (1 - p) * (1 - 6 * p + 6 * p^2)
  }

  list(l1 = l1, l2 = l2, l3 = l3, l4 = l4, p = p)
}

# family constructor

#' Zero-inflated (mixture) Poisson location-scale model family
#'
#' @description
#' The \code{ziPlss} family implements a zero-inflated Poisson model in which one linear predictor controls the
#' Poisson mean and the other controls the probability that an observation is a structural zero. Useable only
#' with \code{\link[mgcv]{gam}}, the linear predictors are specified via a list of formulae. It should be used with
#' care: simply having a large number of zeroes is not an indication of zero inflation.
#'
#' \code{ziPll} evaluates the log-likelihood and its derivatives with respect to the two linear predictors. It is
#' used by \code{ziPlss} and is exported only to allow external use; it is not usually called directly.
#'
#' Requires non-negative integer count data.
#'
#' @details
#' \bold{Specifying the model.} \code{gam} is called with a list containing two formulae. The first has the
#' response on the left hand side and the structure of the linear predictor for the Poisson mean on the right hand
#' side. The second is one sided and specifies the linear predictor for the probability of a structural zero.
#'
#' \bold{Parameterisation and likelihood.} With \eqn{\lambda=\exp(\eta_1)}{lambda = exp(eta_1)} and
#' \eqn{p=1/(1+\exp(-\eta_2))}{p = 1 / (1 + exp(-eta_2))},
#' \deqn{P(Y=0)=p+(1-p)e^{-\lambda},\qquad P(Y=y)=(1-p)\,\frac{e^{-\lambda}\lambda^{y}}{y!},\quad y=1,2,\ldots}{P(Y = 0) = p + (1 - p) exp(-lambda), P(Y = y) = (1 - p) exp(-lambda) lambda^y / y!, y = 1, 2, ...}
#' so that \eqn{E(Y)=(1-p)\lambda}{E(Y) = (1 - p) lambda} and
#' \eqn{\mathrm{Var}(Y)=(1-p)\lambda(1+p\lambda)}{Var(Y) = (1 - p) lambda (1 + p lambda)}. The log-likelihood
#' contribution of observation \eqn{i}{i} is
#' \deqn{\ell_i=\log\left(p_i+(1-p_i)e^{-\lambda_i}\right)\ \ (y_i=0),\qquad \ell_i=\log(1-p_i)+y_i\eta_{1i}-\lambda_i-\log y_i!\ \ (y_i>0).}{l_i = log(p_i + (1 - p_i) exp(-lambda_i)) (y_i = 0); l_i = log(1 - p_i) + y_i eta_1i - lambda_i - log(y_i!) (y_i > 0).}
#' Only identity links are available, because the parameterisation is directly in terms of
#' \eqn{\log\lambda}{log(lambda)} and \eqn{\mathrm{logit}\,p}{logit(p)}. Note the difference from
#' \code{\link[mgcv]{ziplss}}, which is a hurdle model parameterised through the complementary log log of the
#' probability of presence.
#'
#' \bold{Fitted values and prediction.} The fitted values for this family are a two column matrix of the linear
#' predictors \eqn{(\log\lambda,\ \mathrm{logit}\,p)}{(log(lambda), logit(p))}. \code{predict.gam} with
#' \code{type = "link"} returns the same two columns, and \code{type = "response"} returns \eqn{E(Y)}{E(Y)} (with
#' standard errors by the delta method if requested).
#'
#' \bold{Residuals.} \code{residuals} supports \code{"deviance"}, \code{"pearson"}, \code{"working"} and
#' \code{"response"}. Deviance residuals use the saturated log-likelihood, zero for \eqn{y=0}{y = 0} and
#' \eqn{y\log y-y-\log y!}{y log(y) - y - log(y!)} for \eqn{y>0}{y > 0}.
#'
#' \bold{Null deviance.} The null deviance assumes that a single probability of a structural zero and a single
#' Poisson mean are estimated (with the probability profiled out).
#'
#' @section Bias correction of the zero-inflation parameter (\code{correct.bias}):
#' In \code{mgcv}'s general family framework the unpenalised coefficients of every linear predictor are estimated
#' jointly with the penalised coefficients as a penalised-likelihood mode, and the Laplace-approximate REML
#' criterion is used only to choose the smoothing parameters. Marginal-likelihood fitters such as \pkg{glmmTMB}
#' integrate out the smooth coefficients and (under REML) the mean-predictor coefficients, but not the
#' unpenalised coefficients of the zero-inflation predictor. With \code{correct.bias = TRUE} (the default) a
#' one-step correction of the unpenalised coefficients of the zero-inflation predictor towards the
#' marginal-likelihood estimate is applied after fitting, in the \code{postproc} step, using \code{n.newton = 0}.
#' The coefficients before correction are stored in the fitted object as \code{joint_coef} and the full output of
#' \code{\link{bias_correct}} as \code{bias.correction}; \code{Vp}, \code{edf}, \code{sp}, \code{aic} and
#' \code{gcv.ubre} are those of the uncorrected (joint-mode) fit. If the correction fails, the uncorrected fit is
#' returned and the error is stored in \code{bc}. The correction assumes smoothing parameter selection by
#' \code{method = "REML"} (or \code{"ML"}); set \code{correct.bias = FALSE} to obtain the uncorrected fit. There is
#' no overdispersion parameter in this family, so the correction is expected to matter less than for
#' \code{\link{ziNBlss}}.
#'
#' @param link a list of two character strings specifying the links. Currently only \code{"identity"} is
#'   available, as the parameterisation is directly in terms of \eqn{\log\lambda}{log(lambda)} and
#'   \eqn{\mathrm{logit}\,p}{logit(p)}.
#' @param correct.bias logical. If \code{TRUE} (default), apply the one-step marginal-likelihood correction of the
#'   unpenalised coefficients of the zero-inflation linear predictor after fitting. See the section on bias
#'   correction and \code{\link{bias_correct}}.
#' @param y response: a vector of non-negative integer counts.
#' @param eta1,eta2 numeric vectors, of the same length as \code{y}, holding the two linear predictors
#'   \eqn{\log\lambda}{log(lambda)} and \eqn{\mathrm{logit}\,p}{logit(p)}.
#' @param deriv integer between 0 and 4 giving the number of derivatives to compute: 0 returns the
#'   log-likelihood only; 1 additionally the first and second derivatives (and the expected Hessian);
#'   2 additionally the third derivatives; 4 additionally the fourth derivatives.
#'
#' @return
#' For \code{ziPlss}, an object inheriting from class \code{general.family}.
#'
#' For \code{ziPll}, a list with elements
#' \describe{
#'   \item{\code{l}}{vector of log-likelihood contributions.}
#'   \item{\code{l1}}{\eqn{n\times 2}{n x 2} matrix of first derivatives with respect to
#'     \eqn{\eta_1,\eta_2}{eta_1, eta_2}.}
#'   \item{\code{l2}}{\eqn{n\times 3}{n x 3} matrix of second derivatives, columns ordered 11, 12, 22.}
#'   \item{\code{l3}}{\eqn{n\times 4}{n x 4} matrix of third derivatives, columns ordered 111, 112, 122, 222.}
#'   \item{\code{l4}}{\eqn{n\times 5}{n x 5} matrix of fourth derivatives, columns ordered
#'     1111, 1112, 1122, 1222, 2222.}
#'   \item{\code{El2}}{\eqn{n\times 3}{n x 3} matrix holding the expected Hessian in the same layout as \code{l2}.}
#' }
#' Elements that are not requested by \code{deriv} are \code{NULL}. All derivatives are with respect to the linear
#' predictors \eqn{\eta_1,\eta_2}{eta_1, eta_2}.
#'
#' @section WARNINGS:
#' Zero inflated models are often over-used. Having many zeroes in the data does not in itself imply zero
#' inflation. Having too many zeroes \emph{given the model mean} may imply zero inflation.
#'
#' @references
#' Lambert, D. (1992) Zero-inflated Poisson regression, with an application to defects in manufacturing.
#' \emph{Technometrics} 34, 1-14.
#'
#' Wood, S.N., N. Pya and B. Saefken (2016) Smoothing parameter and model selection for general smooth
#' models. \emph{Journal of the American Statistical Association} 111, 1548-1575.
#' \doi{10.1080/01621459.2016.1180986}
#'
#' @seealso \code{\link{ziNBlss}}, \code{\link{hurdNBlss}}, \code{\link{bias_correct}},
#'   \code{\link[mgcv]{ziplss}}, \code{\link[mgcv]{gam}}, \code{\link[mgcv]{general.family}}
#'
#' @examples
#' library(mgcv)
#' set.seed(3)
#' n <- 1000
#' x1 <- runif(n); x2 <- runif(n)
#' lambda <- exp(0.5 + sin(2 * pi * x1))
#' p <- plogis(-1 + 2 * (x2 - 0.5))         # P(structural zero)
#' y <- rpois(n, lambda)
#' y[runif(n) < p] <- 0
#'
#' \donttest{
#' b <- gam(list(y ~ s(x1), ~ s(x2)), family = ziPlss(), method = "REML")
#' b$outer.info
#' summary(b)
#' plot(b, pages = 1)
#' }
#' @rdname ziPlss
#' @export
ziPlss <- function (link = list("identity", "identity"), correct.bias = TRUE)
{
  if (length(link) != 2)
    stop("ziPlss requires 2 links specified as character strings")
  okLinks <- list(c("identity"), c("identity"))
  stats <- list()
  param.names <- c("Poisson mean", "binary probability")
  for (i in 1:2) {
    if (link[[i]] %in% okLinks[[i]])
      stats[[i]] <- make.link(link[[i]])
    else stop(link[[i]], " link not available for ", param.names[i],
              " parameter of ziPlss")
    fam <- structure(list(link = link[[i]], canonical = "none",
                          linkfun = stats[[i]]$linkfun, mu.eta = stats[[i]]$mu.eta),
                     class = "family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }
  residuals <- function(object, type = c("deviance", "pearson", "working", "response")) {
    type <- match.arg(type)
    y    <- object$y
    eta1 <- object$fitted.values[, 1]   # log lambda
    eta2 <- object$fitted.values[, 2]   # logit p, p = P(structural zero)
    lam  <- exp(eta1)
    p    <- plogis(eta2)

    Ey  <- (1 - p) * lam
    rsd <- y - Ey
    if (type == "response") return(rsd)

    if (type == "pearson") {
      Vy <- pmax((1 - p) * lam * (1 + p * lam), .Machine$double.eps)
      return(rsd / sqrt(Vy))
    }
    if (type == "working") return(rsd / pmax(Ey, .Machine$double.eps))  # dE/d eta1 = Ey

    ls <- ifelse(y > 0, y * log(pmax(y, 1)) - y - lgamma(y + 1), 0)
    l  <- ziPll(y, eta1, eta2, deriv = 0)$l
    sign(rsd) * sqrt(pmax(0, 2 * (ls - l)))
  }
  predict <- function(family, se = FALSE, eta = NULL, y = NULL,
                      X = NULL, beta = NULL, off = NULL, Vb = NULL) {
    if (is.null(eta)) {
      if (is.null(off)) off <- list(0, 0)
      off[[3]] <- 0
      for (i in 1:2) if (is.null(off[[i]])) off[[i]] <- 0
      lpi  <- attr(X, "lpi")
      eta1 <- drop(X[, lpi[[1]], drop = FALSE] %*% beta[lpi[[1]]] + off[[1]])
      eta2 <- drop(X[, lpi[[2]], drop = FALSE] %*% beta[lpi[[2]]] + off[[2]])
    } else {
      se <- FALSE
      eta1 <- eta[, 1]; eta2 <- eta[, 2]
    }
    lam <- exp(eta1); p <- plogis(eta2)
    fv  <- list((1 - p) * lam)
    if (!se) return(fv)

    G <- matrix(0, nrow(X), ncol(X))
    G[, lpi[[1]]] <- X[, lpi[[1]], drop = FALSE] * ((1 - p) * lam)
    G[, lpi[[2]]] <- X[, lpi[[2]], drop = FALSE] * (-lam * p * (1 - p))
    fv[[2]] <- sqrt(pmax(0, rowSums((G %*% Vb) * G)))
    names(fv) <- c("fit", "se.fit")
    fv
  }
  rd <- function(mu, wt, scale) {
    ## mu: n x 2 matrix of linear predictors (log lambda, logit p)
    if (NCOL(mu) < 2) stop("rd for ziPlss needs the 2 linear predictors")
    n <- nrow(mu)
    y <- rpois(n, exp(mu[, 1]))
    y[runif(n) < plogis(mu[, 2])] <- 0     # structural zeros
    y
  }
  postproc <- expression({

    # TODO: it is possible that we want the bias correction after the subsequent deviance calculation but I haven't given this much thought
    # Bias correction for unpenalised coefficients in the higher-order predictors
    if (object$family$correct.bias) {
      bc <- try(bias_correct(object, G, n.newton = 0L))
      if(!inherits(bc, "try-error")) {
        object$joint_coef <- bc$coefficients.joint
        names(object$joint_coef) <- G$term.names
        object$coefficients <- bc$coefficients
        object$bias.correction <- bc
        object$linear.predictors[] <- bc$eta
        object$fitted.values[]     <- bc$eta
      } else {
        object$bc <- bc
      }
    }

    ls <- function(y) ifelse(y > 0, y * log(pmax(y, 1)) - y - lgamma(y + 1), 0)

    y_obs <- object$y
    n  <- length(y_obs); n0 <- sum(y_obs == 0)
    ybar <- mean(y_obs)

    prof_obj <- function(el) {            # el = log lambda; p profiled out
      lam <- exp(el)
      p   <- min(max((n0 / n - exp(-lam)) / -expm1(-lam), 1e-10), 1 - 1e-10)
      s   <- sum(ziPll(y_obs, rep(el, n), rep(qlogis(p), n), deriv = 0)$l)
      if (is.finite(s)) -s else 1e10
    }
    lnull <- if (ybar > 0) {
      -optimize(prof_obj, c(log(ybar) - 2, log(max(y_obs) + 1)))$objective
    } else 0

    lsat <- sum(ls(y_obs))
    lhat <- sum(ziPll(y_obs, object$fitted.values[, 1], object$fitted.values[, 2],
                      deriv = 0)$l)
    object$null.deviance <- 2 * (lsat - lnull)
    object$deviance      <- 2 * (lsat - lhat)
  })
  ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL,
                  Hi = NULL, R = NULL, offset = NULL, dH = NULL, db = NULL,
                  deriv = FALSE, nt = 1) {
    gamlss.ncv(X, y, wt, nei, beta, family, llf, H = H, Hi = Hi,
               R = R, offset = offset, dH = dH, db = db, deriv = deriv,
               nt = nt)
  }
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL, D = NULL,
                 eta = NULL, ncv = FALSE, sandwich = FALSE) {
    if (is.null(offset))
      offset <- list(0, 0)
    else offset[[3]] <- 0
    for (i in 1:2) if (is.null(offset[[i]]))
      offset[[i]] <- 0
    jj <- attr(X, "lpi")
    if (is.null(eta)) {
      eta <- X[, jj[[1]], drop = FALSE] %*% coef[jj[[1]]] +
        offset[[1]]
      eta1 <- X[, jj[[2]], drop = FALSE] %*% coef[jj[[2]]] +
        offset[[2]]
    }
    else {
      eta1 <- eta[, 2]
      eta <- eta[, 1]
    }
    lambda <- family$linfo[[1]]$linkinv(eta)
    p <- family$linfo[[2]]$linkinv(eta1)
    zl <- ziPll(y, lambda, p, deriv)
    if (deriv) {
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta), family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(lambda), family$linfo[[2]]$d2link(p))
    }
    g3 <- g4 <- 0
    if (deriv > 1) {
      g3 <- cbind(family$linfo[[1]]$d3link(lambda), family$linfo[[2]]$d3link(p))
    }
    if (deriv > 3) {
      g4 <- cbind(family$linfo[[1]]$d4link(lambda), family$linfo[[2]]$d4link(p))
    }
    if (deriv) {
      i2 <- family$tri$i2
      i3 <- family$tri$i3
      i4 <- family$tri$i4
      de <- gamlss.etamu(zl$l1, zl$l2, zl$l3, zl$l4, ig1,
                         g2, g3, g4, i2, i3, i4, deriv - 1)
      ret <- gamlss.gH(X, jj, de$l1, de$l2, i2, l3 = de$l3,
                       i3 = i3, l4 = de$l4, i4 = i4, d1b = d1b, d2b = d2b,
                       deriv = deriv - 1, fh = fh, D = D, sandwich = sandwich)
      if (ncv) {
        ret$l1 <- de$l1
        ret$l2 = de$l2
        ret$l3 = de$l3
      }
    }
    else ret <- list()
    ret$l <- sum(zl$l)
    ret
  }
  sandwich <- function(y, X, coef, wt, family, offset = NULL) {
    ll(y, X, coef, wt, family, offset = NULL, deriv = 1,
       sandwich = TRUE)$lbb
  }
  initialize <- expression({
    n <- rep(1, nobs)
    if (any(y < 0) || !isTRUE(all.equal(y, round(y)))) {
      stop("Non-negative integer response required with ziPlss")
    }
    if (min(y) == 0 && max(y) == 1) stop("Using ziPlss for binary data makes no sense")
    if (all(y == 0)) stop("ziPlss needs some positive counts")
    if (is.null(start)) {
      use.unscaled <- !is.null(attr(E, "use.unscaled"))
      jj <- attr(x, "lpi")

      start <- (function(y, x, E, jj, use.unscaled, n_em = 3L) {
        ## ---- helper (local, so nothing leaks into mgcv's frame) ----
        pen_fit <- function(xx, ee, yy, ww = NULL) {
          if (!is.null(ww)) {
            ww <- pmax(ww, 1e-8); ww <- ww / mean(ww)
            sw <- sqrt(ww); xx <- xx * sw; yy <- yy * sw
          }
          if (use.unscaled) {
            cf <- qr.coef(qr(rbind(xx, ee)), c(yy, rep(0, nrow(ee))))
            cf[!is.finite(cf)] <- 0
            cf
          } else pen.reg(xx, ee, yy)
        }
        lo_mu <- log(1e-3); hi_mu <- log(max(y) + 1)

        x1 <- x[, jj[[1]], drop = FALSE]; e1 <- E[, jj[[1]], drop = FALSE]
        x2 <- x[, jj[[2]], drop = FALSE]; e2 <- E[, jj[[2]], drop = FALSE]
        nn  <- length(y)
        pos <- y > 0

        ## ---- global moment solution for (lambda, p) ----
        ## E[y] = (1-p) lam,  P(y=0) = p + (1-p) exp(-lam)
        ## eliminating p:  1 - (ybar / lam) * (1 - exp(-lam)) = z0
        ybar <- mean(y); z0 <- mean(!pos)
        lam_g <- ybar
        if (z0 > exp(-ybar) + 1e-8 && z0 < 1) {
          lam_g <- tryCatch(
            uniroot(function(l) 1 - ybar / l * (-expm1(-l)) - z0,
                    c(ybar, ybar / (1 - z0) + 1), tol = 1e-8)$root,
            error = function(e) ybar)
        }
        pbar <- min(max(1 - ybar / lam_g, 0.01), 0.99)

        ## ---- initial values ----
        ## mean: log y on positive rows (zeros carry no usable log-scale information)
        cf1 <- if (sum(pos) >= 5) pen_fit(x1[pos, , drop = FALSE], e1, log(y[pos]))
        else pen_fit(x1, e1, rep(log(lam_g), nn))
        eta1 <- pmin(pmax(drop(x1 %*% cf1), lo_mu), hi_mu)
        mu   <- exp(eta1)

        ## p: constant start from the moment solution
        cf2  <- pen_fit(x2, e2, rep(qlogis(pbar), nn))
        eta2 <- pmin(pmax(drop(x2 %*% cf2), -10), 10)

        ## ---- a few EM-style sweeps (no likelihood evaluations) ----
        for (k in seq_len(n_em)) {
          p  <- plogis(eta2)
          z  <- ifelse(pos, 0, p / (p + (1 - p) * exp(-mu)))   # P(structural | y)

          ## p: IRLS step on soft labels z
          w2   <- pmax(p * (1 - p), 1e-6)
          cf2  <- pen_fit(x2, e2, eta2 + (z - p) / w2, w2)
          eta2 <- pmin(pmax(drop(x2 %*% cf2), -10), 10)

          ## lambda: Poisson IRLS step, rows weighted by P(Poisson-generated)
          cf1  <- pen_fit(x1, e1, eta1 + (y - mu) / mu, (1 - z) * mu)
          eta1 <- pmin(pmax(drop(x1 %*% cf1), lo_mu), hi_mu)
          mu   <- exp(eta1)
        }

        st <- rep(0, ncol(x))
        st[jj[[1]]] <- cf1
        st[jj[[2]]] <- cf2
        st
      })(y, x, E, jj, use.unscaled)
    }
  })
  structure(list(family = "ziPlss", ll = ll, link = paste(link), correct.bias = correct.bias,
                 nlp = 2, ncv = ncv, sandwich = sandwich, tri = trind.generator(2),
                 initialize = initialize, postproc = postproc, residuals = residuals,
                 rd = rd, predict = predict, linfo = stats, d2link = 1,
                 d3link = 1, d4link = 1, ls = 1, available.derivs = 2),
            class = c("general.family", "extended.family", "family"))
}

################################################################################
#### Post-processing bias correction  ##########################################
################################################################################

#' One-step marginal-likelihood (REML-type) correction of unpenalised coefficients
#'
#' @description
#' Moves the unpenalised coefficients of selected linear predictors of a \code{\link{hurdNBlss}},
#' \code{\link{ziNBlss}} or \code{\link{ziPlss}} fit from their joint penalised-likelihood mode (as computed by
#' \code{\link[mgcv]{gam}}) towards the estimate that maximises a marginal (Laplace/REML-type) likelihood, in which
#' the smooth and mean-predictor coefficients are integrated out but the unpenalised coefficients of the
#' zero-probability and dispersion predictors are not. This is the convention of \pkg{glmmTMB} and \pkg{sdmTMB}
#' with \code{REML = TRUE}. The correction is a single Newton-type step computed from the fitted model, with no
#' refitting, and is run automatically after fitting when a family is called with \code{correct.bias = TRUE}.
#'
#' @details
#' \bold{Why a correction is needed.} For a general family \code{gam} estimates all coefficients, including
#' unpenalised ones such as the intercept of the overdispersion predictor, jointly by maximising the penalised
#' log-likelihood \eqn{\ell_p}{l_p} for given smoothing parameters; the Laplace-approximate REML criterion
#' only chooses the smoothing parameters. A marginal-likelihood fit instead maximises, over the unpenalised
#' parameters \eqn{\psi}{psi}, the criterion
#' \deqn{Q(\psi)=\ell_p^{*}(\psi)-\tfrac12\log\left|\mathbf H_{RR}(\hat R_\psi,\psi)\right|,}{Q(psi) = l_p*(psi) - 0.5 log|H_RR(R_hat(psi), psi)|,}
#' where \eqn{R}{R} denotes the integrated-out coefficients, \eqn{\hat R_\psi}{R_hat(psi)} their penalised mode given
#' \eqn{\psi}{psi}, \eqn{\ell_p^{*}(\psi)=\ell_p(\hat R_\psi,\psi)}{l_p*(psi) = l_p(R_hat(psi), psi)} the profile, and
#' \eqn{\mathbf H_{RR}=-\partial^2\ell_p/\partial R\,\partial R^{\top}}{H_RR = -d^2 l_p / dR dR'} the penalised
#' Hessian for \eqn{R}{R}. The joint mode maximises the profile \eqn{\ell_p^{*}}{l_p*} alone; the log-determinant
#' (a Cox-Reid type adjustment, analogous to the \eqn{n-p}{n - p} divisor in REML estimation of a Gaussian variance)
#' is what separates the two. Because the information about the mean increases with the dispersion parameter, the
#' joint mode overestimates it, increasingly so as the effective degrees of freedom of the mean predictor grow.
#'
#' \bold{Which coefficients.} For each linear predictor to be corrected, the coefficients are rotated into the
#' null space of its penalty (\eqn{\psi}{psi}: for example the intercept and linear terms of a smooth, or all the
#' coefficients of an unpenalised predictor) and its penalised directions (which belong to \eqn{R}{R}). The mean
#' predictor is always integrated out and is never corrected. The default sets follow the \pkg{glmmTMB} REML
#' convention: for \code{ziNBlss} the zero-inflation and dispersion predictors (\code{c(FALSE, TRUE, TRUE)}); for
#' \code{hurdNBlss} the dispersion predictor only (\code{c(FALSE, FALSE, TRUE)}; the zero-probability predictor
#' separates from the rest and gives the same dispersion shift whether or not it is included); for \code{ziPlss}
#' the zero-inflation predictor (\code{c(FALSE, TRUE)}).
#'
#' \bold{The one-step shift.} A Newton step on \eqn{Q}{Q} from the joint mode, where the score of \eqn{\ell_p^{*}}{l_p*}
#' is zero, gives
#' \deqn{\tilde\psi=\hat\psi-\tfrac12\,\mathbf V_{\psi\psi}\,\mathbf D,\qquad D_k=\frac{d}{d\psi_k}\log\left|\mathbf H_{RR}(\hat R_\psi,\psi)\right|,}{psi_tilde = psi_hat - 0.5 V_psipsi D, D_k = d log|H_RR(R_hat(psi), psi)| / d psi_k,}
#' where \eqn{\mathbf V_{\psi\psi}}{V_psipsi} is the \eqn{\psi\psi}{psi psi} block of the inverse of the full penalised
#' Hessian (the \code{Vp} matrix of the fit, rotated), i.e. the inverse of the information about \eqn{\psi}{psi} after
#' profiling out \eqn{R}{R}, and the curvature of the log-determinant is neglected. With
#' \eqn{\ell_{i,abc}}{l_i,abc} the third derivatives of the log-likelihood contributions with respect to the linear
#' predictors, \eqn{h_{ab,i}=\mathbf a_{a,i}^{\top}\mathbf H_{RR}^{-1}\mathbf a_{b,i}}{h_ab,i = a_a,i' H_RR^-1 a_b,i}
#' (with \eqn{\mathbf a_{c,i}}{a_c,i} the row of the design of predictor \eqn{c}{c} in the coordinates of \eqn{R}{R}),
#' and \eqn{e^{(k)}_{c,i}=d\eta_{ci}/d\psi_k}{e_c,i = d eta_ci / d psi_k},
#' \deqn{D_k=-\sum_i\sum_{a,b}\sum_c \ell_{i,abc}\,e^{(k)}_{c,i}\,h_{ab,i}.}{D_k = - sum_i sum_a,b sum_c l_i,abc e_c,i h_ab,i.}
#' The term \eqn{e^{(k)}_{c,i}}{e_c,i} has a direct part (the design of \eqn{\psi_k}{psi_k} in predictor \eqn{c}{c})
#' and an indirect part through the dependence of \eqn{\hat R_\psi}{R_hat(psi)} on \eqn{\psi}{psi}, which comes from
#' differentiating the penalised score. These give the components \code{D_direct} and \code{D_indirect}. For a scalar
#' dispersion parameter the direct part is proportional to the effective degrees of freedom of the mean predictor,
#' which explains why the joint-mode bias grows with the flexibility of the mean surface.
#'
#' \bold{Refining the integrated coefficients.} By default only \eqn{\psi}{psi} is moved and \eqn{R}{R} is left at its
#' joint-mode value, so the corrected coefficient vector is not exactly a stationary point of \eqn{\ell_p}{l_p}
#' with respect to \eqn{R}{R}. Setting \code{first.order.R = TRUE} also moves \eqn{R}{R} by its first-order response
#' to the shift, and \code{n.newton > 0} then re-solves \eqn{R}{R} at the corrected \eqn{\psi}{psi} by chord Newton
#' iterations that re-use the Hessian from the joint solution, so each iteration costs one gradient evaluation.
#' \code{n.newton} is a maximum: iterations stop when the gradient norm falls below \code{newton.tol}. Smoothing
#' parameters are not re-optimised.
#'
#' \bold{Diagnostics.} The first-order step is accurate if the shift is small relative to the curvature of
#' \eqn{Q}{Q}. \code{step_sd} reports the size of the shift in units of standard deviations of the profile
#' likelihood, \eqn{\tfrac12\sqrt{\mathbf D^{\top}\mathbf V_{\psi\psi}\mathbf D}}{0.5 sqrt(D' V_psipsi D)}, and a warning is
#' issued if it exceeds \code{step.warn}. \code{pred_gain} is the predicted increase in \eqn{Q}{Q}; with
#' \code{criterion = TRUE} the realised increase is also computed (meaningful once \eqn{R}{R} has been re-solved, e.g.
#' with \code{n.newton >= 3}) and \code{gain_ratio} is the ratio of realised to predicted gain, which is close to
#' one when the quadratic approximation is adequate. \code{cross_share} gives, for each shifted coordinate, the
#' fraction of its shift that comes through the other shifted coordinates (large when, for example, the
#' zero-inflation and dispersion parameters are strongly correlated).
#'
#' \bold{Use inside \code{postproc}.} When called from a family's \code{postproc} step, \code{m} is the partially
#' assembled fitted object and \code{G} the model set-up object, from which the design matrix (with its
#' \code{"lpi"} attribute), the smooth definitions, the model frame, the offsets and the term names are taken.
#' Offsets enter through the linear predictors. The call made by the families is
#' \code{bias_correct(object, G, n.newton = 0L)}.
#'
#' @section Limitations:
#' The derivation assumes smoothing parameters selected by REML or ML, identity links for all linear predictors,
#' and, as implemented, a single set of coefficients per predictor with no dropped (\code{bdrop}) coefficients;
#' penalties supplied through \code{paraPen} are not handled. The design matrices in the coordinates of the
#' integrated and unintegrated coefficients are stored densely, so memory grows like \eqn{n}{n} times the number of
#' coefficients times the number of linear predictors. The shift is first order and may under- or over-correct when it
#' is large (see \code{step_sd} and \code{gain_ratio}). The correction is available for the hurdle and mixture
#' negative binomial and the zero-inflated Poisson families of this package only. The estimated covariance matrix
#' \code{Vp} refers to the joint-mode fit; intervals for corrected coefficients are formed by re-centring on the
#' corrected estimate.
#'
#' @param m a fitted \code{\link[mgcv]{gam}} object from \code{\link{hurdNBlss}}, \code{\link{ziNBlss}} or
#'   \code{\link{ziPlss}} (or, inside \code{postproc}, the partially assembled object). It must contain the
#'   coefficients, \code{Vp}, \code{sp}, the smooth definitions, the response \code{y} and the family.
#' @param G \code{NULL} (default), in which case the design matrix is obtained from
#'   \code{predict(m, type = "lpmatrix")}, or the \code{gam} set-up object used inside \code{postproc}
#'   (containing \code{X}, \code{smooth}, \code{mf}, \code{offset} and \code{term.names}).
#' @param apply_to logical vector with one entry per linear predictor indicating which of the computed shifts to apply.
#'   The default applies all shifts that were computed (the family default set, see Details). Shifts can only be
#'   applied to predictors that were corrected. The \emph{computed} shifts do not depend on this argument; it only
#'   selects which are added to the coefficients, which is useful for diagnostics.
#' @param n.newton maximum number of chord Newton iterations used to re-solve the integrated coefficients at the
#'   corrected parameters. \code{0} leaves them at their joint-mode values.
#' @param first.order.R logical. If \code{TRUE}, move the integrated coefficients by their first-order response to the
#'   shift before any Newton iterations.
#' @param criterion logical. If \code{TRUE}, evaluate the marginal criterion \eqn{Q}{Q} (and its components) at the
#'   joint and corrected coefficients and compute \code{gain_ratio}.
#' @param tol relative tolerance used to identify the null space of a penalty (eigenvalues no larger than
#'   \code{tol} times the largest eigenvalue are treated as zero).
#' @param step.warn warn if the shift exceeds this many standard deviations of the profile likelihood
#'   (\code{step_sd}).
#' @param diagnostics logical. If \code{TRUE} (default), compute the gradient diagnostics (\code{grad_R}) even when
#'   no refinement is requested.
#' @param newton.tol gradient norm at which the Newton iterations stop. The default is ten times the gradient norm
#'   at the joint solution (the numerical floor), with a lower bound of \code{1e-10}.
#' @param verbose logical. If \code{TRUE}, report which predictors are corrected and which shifts are applied.
#'
#' @return
#' A list with components
#' \describe{
#'   \item{\code{coefficients}}{the corrected coefficient vector.}
#'   \item{\code{coefficients.joint}}{the original joint-mode coefficient vector.}
#'   \item{\code{delta}}{\code{coefficients - coefficients.joint}, in the original coordinates (the shift lies in the
#'     null space of the penalty, so for a smooth it may involve the last, unpenalised, coefficients as well as the
#'     intercept).}
#'   \item{\code{shift}, \code{shift_applied}, \code{apply_to}}{the computed shift of the unpenalised coordinates
#'     (in rotated coordinates, whose orientation within a null space is arbitrary), the shift actually applied, and
#'     the logical vector selecting predictors.}
#'   \item{\code{shift_by_predictor}}{\code{shift} split by linear predictor.}
#'   \item{\code{D_direct}, \code{D_indirect}}{the direct and indirect parts of the derivative \eqn{D}{D} of the
#'     log-determinant.}
#'   \item{\code{V_psi}}{the covariance matrix \eqn{\mathbf V_{\psi\psi}}{V_psipsi} of the unintegrated coordinates.}
#'   \item{\code{step_sd}, \code{pred_gain}, \code{gain_ratio}, \code{cross_share}}{diagnostics for the first-order
#'     step (see Details); \code{gain_ratio} is \code{NA} unless \code{criterion = TRUE}.}
#'   \item{\code{nN}, \code{nU}}{the number of unintegrated and integrated coordinates contributed by each linear
#'     predictor.}
#'   \item{\code{eta}}{matrix of the corrected linear predictors (including offsets), one column per predictor.}
#'   \item{\code{grad_R}, \code{grad_trace}, \code{newton_rate}, \code{newton_iter}}{gradient diagnostics for the
#'     integrated coefficients: \code{grad_R} holds the gradient norm at the joint solution (\code{joint}), after the
#'     shift and any first-order update (\code{shifted}) and after the Newton iterations (\code{final});
#'     \code{grad_trace} the norm at the start and after each iteration; \code{newton_rate} the geometric-mean
#'     contraction per iteration; \code{newton_iter} the number of iterations used.}
#'   \item{\code{criterion}}{if \code{criterion = TRUE}, a matrix with rows \code{joint} and \code{corrected} and
#'     columns \code{ll_pen} (penalised log-likelihood), \code{logdet_HRR} and \code{Q}.}
#'   \item{\code{edf_mu}}{effective degrees of freedom of the mean predictor.}
#' }
#' If there are no unpenalised coefficients to correct, only \code{coefficients}, \code{coefficients.joint} and
#' \code{shift} (of length zero) are returned.
#'
#' @references
#' Cox, D.R. and N. Reid (1987) Parameter orthogonality and approximate conditional inference. \emph{Journal of the
#' Royal Statistical Society B} 49, 1-39.
#'
#' Wood, S.N. (2011) Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric
#' generalized linear models. \emph{Journal of the Royal Statistical Society B} 73, 3-36.
#'
#' Wood, S.N., N. Pya and B. Saefken (2016) Smoothing parameter and model selection for general smooth models.
#' \emph{Journal of the American Statistical Association} 111, 1548-1575. \doi{10.1080/01621459.2016.1180986}
#'
#' Brooks, M.E. et al. (2017) glmmTMB balances speed and flexibility among packages for zero-inflated generalized
#' linear mixed modeling. \emph{The R Journal} 9, 378-400.
#'
#' @seealso \code{\link{hurdNBlss}}, \code{\link{ziNBlss}}, \code{\link{ziPlss}}, \code{\link[mgcv]{gam}}
#'
#' @examples
#' library(mgcv)
#' set.seed(2)
#' n <- 1000
#' x1 <- runif(n); x2 <- runif(n)
#' mu  <- exp(0.5 + sin(2 * pi * x1))
#' p   <- plogis(-1 + 2 * (x2 - 0.5))
#' y <- rnbinom(n, mu = mu, size = 2)
#' y[runif(n) < p] <- 0
#'
#' \donttest{
#' ## fit without the automatic correction, then apply it explicitly
#' b0 <- gam(list(y ~ s(x1), ~ s(x2), ~ 1), family = ziNBlss(correct.bias = FALSE),
#'           method = "REML")
#' bc <- bias_correct(b0)                       # shift of the unpenalised coefficients only
#' bc$delta[bc$delta != 0]                      # which coefficients moved
#' bc$step_sd                                   # size of the step in SD units
#'
#' ## also re-solve the integrated coefficients (chord Newton) and evaluate the criterion
#' bc2 <- bias_correct(b0, n.newton = 20, criterion = TRUE)
#' bc2$grad_R; bc2$newton_iter; bc2$criterion
#' }
#' @export
bias_correct <- function(m, G = NULL, apply_to = NULL,
                         n.newton = 20L, first.order.R = FALSE,
                         criterion = FALSE, tol = 1e-8, step.warn = 2, diagnostics = TRUE,
                         newton.tol = NULL, verbose = FALSE) {

  ## defaults for which linear predictors are un-integrated (glmmTMB-REML rule, fixed per family)
  sets <- list(ziNBlss   = c(FALSE, TRUE, TRUE),
               hurdNBlss = c(FALSE, FALSE, TRUE),
               ziPlss    = c(FALSE, TRUE))
  lls  <- c(ziNBlss = "ziNBll", hurdNBlss = "hurdNBll", ziPlss = "ziPll")

  fam <- m$family$family
  if (!fam %in% names(sets)) stop("bias correction is not defined for family '", fam, "'")
  ll_fun  <- get(lls[[fam]], mode = "function")
  y       <- m$y
  correct <- sets[[fam]]

  if (is.null(apply_to)) {
    apply_to <- correct
  } else {
    if (length(apply_to) != length(correct))
      stop("incorrect number of predictors specified in apply_to for the model's family")
    if (!is.logical(apply_to))
      stop("apply_to must be a logical vector indicating the linear predictors whose shifts are applied")
  }
  apply_to <- apply_to & correct            # can only apply shifts that were computed
  if (verbose) message("correct: ", paste(correct, collapse = " "),
                       "; apply_to: ", paste(apply_to, collapse = " "))

  if (any(m$bdrop)) stop("coefficients were dropped in the fit (bdrop); bias correction not supported")

  ## model matrix and coefficient blocks

  if (is.null(G)) {
    Xp <- predict(m, type = "lpmatrix")   # outside postproc
    lp  <- attr(Xp, "lpi")
    if (is.null(lp)) stop("the model matrix must carry the \"lpi\" attribute (e.g, use G$X)")
    nlp <- length(lp)
  } else {
    Xp <- G$X
    lp  <- attr(Xp, "lpi")
    if (is.null(lp)) stop("the model matrix must carry the \"lpi\" attribute (e.g., use G$X)")
    nlp <- length(lp)

    for (k in seq_along(G$smooth)) {
      sm <- G$smooth[[k]]

      ii <- sm$first.para:sm$last.para

      ## Rebuild the smooth basis as predict.gam does
      Xp[, ii] <- mgcv::PredictMat(sm, G$mf)
    }

    ## Initialise one offset per linear predictor
    tmp.offset <- rep(list(0), nlp)

    if (!all(sapply(G$offset, is.null))) {
      for (ii in 1:length(G$offset)) {
        if (!is.null(G$offset[[ii]])) {
          tmp.offset[[ii]] <- G$offset[[ii]]
        }
      }
    }
    attr(Xp, "model.offset") <- tmp.offset
    ## Match the lpmatrix attributes
    colnames(Xp) <- G$term.names
  }

  off <- attr(Xp, "model.offset")
  n   <- nrow(Xp)
  if (length(correct) != nlp) stop("'correct' must have one entry per linear predictor")
  betas <- m$coefficients
  X   <- lapply(lp, function(i) Xp[, i, drop = FALSE])
  q   <- lengths(lp)
  idx <- unlist(lp)                           # coefficient order used throughout

  ## helper: block-diagonal assembly of (nr[j] x nc[j]) blocks
  blk <- function(L, nr, nc) {
    M <- matrix(0, sum(nr), sum(nc)); r0 <- 0; c0 <- 0
    for (j in seq_along(L)) {
      if (nr[j] > 0 && nc[j] > 0) M[r0 + seq_len(nr[j]), c0 + seq_len(nc[j])] <- L[[j]]
      r0 <- r0 + nr[j]; c0 <- c0 + nc[j]
    }
    M
  }

  ## per-predictor penalty blocks: sp-free (for the null-space split) and sp-weighted (penalty).
  ## The smoothing parameter of each penalty is located through first.sp (robust to te()/multi-penalty terms)
  S_unit <- S_sp <- lapply(q, function(qj) matrix(0, qj, qj))
  for (sm in m$smooth) {
    ci <- sm$first.para:sm$last.para
    j  <- which(vapply(lp, function(l) all(ci %in% l), logical(1)))[1]
    for (l in seq_along(sm$S)) {
      ii <- match(ci, lp[[j]])
      S_unit[[j]][ii, ii] <- S_unit[[j]][ii, ii] + sm$S[[l]]
      S_sp[[j]][ii, ii]   <- S_sp[[j]][ii, ii] + m$sp[sm$first.sp + l - 1] * sm$S[[l]]
    }
  }

  ## split each corrected predictor into penalised (R) and unpenalised (psi) directions
  U <- N <- vector("list", nlp)
  for (j in seq_len(nlp)) {
    if (correct[j]) {
      ev  <- eigen(S_unit[[j]], symmetric = TRUE)
      is0 <- ev$values <= tol * max(ev$values, 1e-300)
      N[[j]] <- ev$vectors[,  is0, drop = FALSE]
      U[[j]] <- ev$vectors[, !is0, drop = FALSE]
    } else {
      U[[j]] <- diag(q[j]); N[[j]] <- matrix(0, q[j], 0)
    }
  }
  nU <- vapply(U, ncol, 1L); nN <- vapply(N, ncol, 1L)
  nR <- sum(nU); m_ps <- sum(nN)
  if (m_ps == 0) {
    message("nothing to correct (no unpenalised coefficients in the selected predictors)")
    return(list(coefficients = betas, coefficients.joint = betas, shift = numeric(0)))
  }
  Tm  <- cbind(blk(U, q, nU), blk(N, q, nN))   # gamma_cat = Tm %*% c(u_R, psi), orthogonal
  iR  <- seq_len(nR); iP <- nR + seq_len(m_ps)

  ## designs of each predictor in R and psi coordinates (n x nR, n x m_ps)
  offR <- c(0, cumsum(nU)); offP <- c(0, cumsum(nN))
  A <- Z <- vector("list", nlp)
  for (j in seq_len(nlp)) {
    A[[j]] <- matrix(0, n, nR); Z[[j]] <- matrix(0, n, m_ps)
    if (nU[j] > 0) A[[j]][, offR[j] + seq_len(nU[j])] <- X[[j]] %*% U[[j]]
    if (nN[j] > 0) Z[[j]][, offP[j] + seq_len(nN[j])] <- X[[j]] %*% N[[j]]
  }

  ## penalised Hessian in (R, psi) coordinates; V_psipsi via the Schur complement
  Hf   <- chol2inv(chol(m$Vp))[idx, idx, drop = FALSE]
  Hn   <- crossprod(Tm, Hf %*% Tm)
  Hi   <- chol2inv(chol(Hn[iR, iR, drop = FALSE]))                 # H_RR^{-1}
  Spsi <- Hn[iP, iP, drop = FALSE] - Hn[iP, iR, drop = FALSE] %*% Hi %*% Hn[iR, iP, drop = FALSE]
  Vpsi <- solve(Spsi)

  get_eta <- function(gam) lapply(seq_len(nlp), function(j) drop(X[[j]] %*% gam[lp[[j]]]) + off[[j]])
  call_ll <- function(eta, deriv) do.call(ll_fun, c(list(y), eta, list(deriv = deriv)))
  tri <- mgcv:::trind.generator(nlp)
  zl  <- call_ll(get_eta(betas), 2)
  l2f <- function(ia, ib)     zl$l2[, tri$i2[ia, ib]]
  l3f <- function(ia, ib, ic) zl$l3[, tri$i3[ia, ib, ic]]
  nz3 <- colSums(abs(zl$l3)) > 0                                   # skip identically-zero derivatives
  hasR <- which(nU > 0)

  HA <- lapply(seq_len(nlp), function(ia) A[[ia]] %*% Hi)
  Hab <- list()
  for (ia in hasR) for (ib in hasR) Hab[[paste(ia, ib)]] <- rowSums(HA[[ia]] * A[[ib]])

  Dfun <- function(e) {                                            # -sum l_abc e_c h_ab
    tot <- 0
    for (ia in hasR) for (ib in hasR) for (ic in seq_len(nlp)) {
      if (!nz3[tri$i3[ia, ib, ic]]) next
      tot <- tot + sum(l3f(ia, ib, ic) * e[[ic]] * Hab[[paste(ia, ib)]])
    }
    -tot
  }

  D_dir <- D_ind <- numeric(m_ps); dRmat <- matrix(0, nR, m_ps)
  for (kk in seq_len(m_ps)) {
    zk <- lapply(seq_len(nlp), function(j) Z[[j]][, kk])
    g  <- numeric(nR)
    for (ia in hasR) {
      sc <- 0
      for (ic in seq_len(nlp)) sc <- sc + l2f(ia, ic) * zk[[ic]]
      g <- g + drop(crossprod(A[[ia]], sc))
    }
    dR <- drop(Hi %*% g); dRmat[, kk] <- dR                         # dR_hat / d psi_k
    eR <- lapply(seq_len(nlp), function(j) drop(A[[j]] %*% dR))     # indirect d eta_c / d psi_k
    D_dir[kk] <- Dfun(zk)
    D_ind[kk] <- Dfun(eR)
  }
  D    <- D_dir + D_ind
  VD   <- drop(Vpsi %*% D)
  dpsi <- 0.5 * VD                                                  # one-step shift of psi

  ## diagnostics for the first-order step
  step_sd     <- sqrt(0.25 * sum(D * VD))
  pred_gain   <- 0.125 * sum(D * VD)
  cross_share <- (VD - diag(Vpsi) * D) / VD
  if (step_sd > step.warn)
    warning(sprintf("one-step shift is %.2f SD in the information metric (> %.1f): the quadratic approximation may be unreliable",
                    step_sd, step.warn))

  ## corrected coordinates
  theta  <- drop(crossprod(Tm, betas[idx]))
  uR     <- theta[iR]; psi <- theta[iP]
  mask   <- apply_to[rep(seq_len(nlp), nN)]                         # which psi coordinates are updated
  dpsi_a <- dpsi * mask                                             # shift actually applied
  psi_n  <- psi - dpsi_a
  uR_n   <- if (first.order.R) uR - drop(dRmat %*% dpsi_a) else uR

  Sbd <- blk(S_sp, q, q)                                            # sp-weighted penalty, coefficient order
  score_R <- function(uRv, psiv) {                                  # penalised score for R at (uRv, psiv)
    gcat  <- drop(Tm %*% c(uRv, psiv))
    gfull <- betas; gfull[idx] <- gcat
    z1 <- call_ll(get_eta(gfull), 1)
    sc <- unlist(lapply(seq_len(nlp), function(j) drop(crossprod(X[[j]], z1$l1[, j]))))
    drop(crossprod(Tm, sc - drop(Sbd %*% gcat)))[iR]
  }

  ## gradient diagnostics and chord-Newton polish (n.newton is a maximum)
  grad <- c(joint = NA_real_, shifted = NA_real_, final = NA_real_)
  trace <- NA_real_; rate <- NA_real_; n_done <- 0L
  if (diagnostics || first.order.R || n.newton > 0) {
    grad["joint"] <- sqrt(sum(score_R(uR, psi)^2))                  # ~0 at the joint solution (the floor)
    scv   <- score_R(uR_n, psi_n)
    trace <- sqrt(sum(scv^2))
    tol_n <- if (is.null(newton.tol)) max(10 * grad[["joint"]], 1e-10) else newton.tol
    for (it in seq_len(n.newton)) {
      if (trace[length(trace)] <= tol_n) break
      uR_n  <- uR_n + drop(Hi %*% scv)                              # chord Newton, fixed H_RR
      scv   <- score_R(uR_n, psi_n)
      trace <- c(trace, sqrt(sum(scv^2)))
    }
    n_done <- length(trace) - 1
    grad["shifted"] <- trace[1]
    grad["final"]   <- trace[length(trace)]
    if (n_done > 0 && trace[1] > 0 && trace[length(trace)] > 0)
      rate <- (trace[length(trace)] / trace[1])^(1 / n_done)
  }

  gcat_n <- drop(Tm %*% c(uR_n, psi_n))
  coef_n <- betas; coef_n[idx] <- gcat_n
  eta_n  <- do.call(cbind, get_eta(coef_n))

  ## glmmTMB-type marginal criterion at a given coefficient vector (R should be at its conditional mode)
  crit_at <- function(coefs) {
    z2   <- call_ll(get_eta(coefs), 2)
    gcat <- coefs[idx]
    pen  <- 0.5 * drop(crossprod(gcat, Sbd %*% gcat))
    HRR  <- crossprod(Tm, Sbd %*% Tm)[iR, iR, drop = FALSE]
    for (ia in hasR) for (ib in hasR)
      HRR <- HRR - crossprod(A[[ia]], A[[ib]] * z2$l2[, tri$i2[ia, ib]])
    ld   <- as.numeric(determinant(HRR, logarithm = TRUE)$modulus)
    c(ll_pen = sum(z2$l) - pen, logdet_HRR = ld, Q = sum(z2$l) - pen - 0.5 * ld)
  }
  crit <- NULL; gain_ratio <- NA_real_
  if (criterion) {
    crit <- rbind(joint = crit_at(betas), corrected = crit_at(coef_n))
    gain_ratio <- (crit["corrected", "Q"] - crit["joint", "Q"]) / pred_gain
  }

  list(coefficients = coef_n, coefficients.joint = betas,
       delta = coef_n - betas,
       shift = dpsi, shift_applied = dpsi_a, apply_to = apply_to,
       shift_by_predictor = split(dpsi, rep(seq_len(nlp), nN)),
       D_direct = D_dir, D_indirect = D_ind, V_psi = Vpsi,
       step_sd = step_sd, pred_gain = pred_gain, gain_ratio = gain_ratio, cross_share = cross_share,
       nN = nN, nU = nU, eta = eta_n, grad_R = grad, grad_trace = trace, newton_rate = rate,
       newton_iter = n_done, criterion = crit,
       edf_mu = sum(m$edf[lp[[1]]]))
}
