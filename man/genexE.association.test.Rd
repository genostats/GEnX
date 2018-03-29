\name{genexE.association.test}
\alias{genexE.association.test}

\title{ Gene x Environement Association Test }

\description{ Association tests between phenotype and interaction between SNPs and environment variable. }
\usage{ 
genexE.association.test(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), E,
                 method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                 test = c("score", "wald", "lrt"), df=c(1,2,3), K, eigenK, beg = 1, 
                 end = ncol(x), p = 0, tol = .Machine$double.eps^0.25, ...)
}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{Y}{ The phenotype vector. Default is the column (\code{pheno}) of \code{x@ped} }
  \item{X}{ A covariable matrix. The default is a column vector of ones, to include an intercept in the model }
  \item{E}{ The environmental variable vector }
  \item{method}{ Method to use: \code{"lm"} for (generalized) linear model, and \code{"lmm"} for (generalized) linear mixed model }
  \item{response}{ Is \code{"Y"} a quantitative or a binary phenotype?}
  \item{test}{ Which test to use. For binary phenotypes, \code{test = "score"} is mandatory}
  \item{df}{ Which null hypothesis to use. See details}
  \item{K}{ A Genetic Relationship Matrix (as produced by \code{\link{GRM}}), or a list of such matrices. For \code{test = "score"}. }
  \item{eigenK}{ Eigen decomposition of the Genetic Relationship Matrix (as produced by the function \code{eigen}). 
                 For \code{test = "wald"} or \code{"lrt"}. }
  \item{beg}{ Index of the first SNP tested for association }
  \item{end}{ Index of the last SNP tested for association }
  \item{p}{ Number of Principal Components to include in the model with fixed effect (for \code{test = "wald"} or \code{"lrt"}) }
  \item{tol}{ Parameter for the likelihood maximization (as in \code{optimize}) }
  \item{...}{ Additional parameters for \code{\link{lmm.aireml}} or \code{\link{logistic.mm.aireml}} (if \code{test = "score"}). }
}

\details{
  Tests the association between the phenotype and interaction between environmental variable and requested SNPs in \code{x}. At most,
  3 possible null hypothesis are availaible (depending on test and response) :
  - for \code{df=1},
  \deqn{ H_0 : beta_{ExSNP}=0 }{ H_0 : beta_ExSNP=0 },
  - for \code{df=2},
  \deqn{ H_0 : beta_{SNP}=beta_{ExSNP}=0 }{ H_0 : beta_SNP=beta_ExSNP=0 },
  - for \code{df=3},
  \deqn{ H_0 : beta_E=beta_{SNP}=beta_{ExSNP}=0 }{ H_0 : beta_E=beta_SNP=beta_ExSNP=0 }.  

  If \code{method = "lm"} and \code{response = "quantitative"} are used, a simple linear regression 
  is performed to test each SNP in the considered interval. Precisely, the following model is
  considered for each SNP,
  \deqn{ Y = (X|PC)\alpha + E\beta_E + G\beta_{SNP} + E*SNP\beta_{ExSNP} + \varepsilon }{ Y = (X|PC) alpha + E beta_E + SNP beta_SNP + E*SNP beta_ExSNP + epsilon }
  with \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)},
  \eqn{E} the vector of the environmental variable, 
  \eqn{SNP} the genotype vector of the SNP, 
  \eqn{X} the covariates matrix, and \eqn{PC} the matrix of the first \eqn{p} principal components.
  A Wald test for given null hypothesis is performed (independently of the value of \code{test}).

  If\code{method = "lm"} and \code{response = "binary"}, a similar model is used for a logistic 
  regression (Wald test).
 
  If \code{method = "lmm"} and \code{response = "quantitative"}, the following model in considered for each SNP
  \deqn{ Y = (X|PC)\alpha + E\beta_E + G\beta_{SNP} + E*SNP\beta_{ExSNP} + \omega + \varepsilon }{ Y = (X|PC) alpha + E beta_E + SNP beta_SNP + E*SNP beta_ExSNP + omega + epsilon }
  with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)} and \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
  \eqn{E} is the vector of the environmental variable, \eqn{SNP} is the genotype vector of the SNP, \eqn{K} is a Genetic Relationship Matrix (GRM)
  \eqn{X} the covariates matrix, and \eqn{PC} the matrix of the first \eqn{p} principal components. The 3 possible null hypothysis are implemented.

  If \code{test = "score"}, all parameters are estimated with the same procedure as in 
  \code{\link{lmm.aireml}} and the argument \code{K} is used to specify the GRM matrix (or a list of GRM
  matrices for inclusion of several random effects in the model). If \code{p} is positive, the paramater \code{eigenK}
  needs to be given as well.
  For Wald and LRT tests the procedure used is the same as in \code{\link{lmm.diago}} and \code{eigenK} is used to 
  specify the GRM matrix.

  If \code{method = "lmm"} and \code{response = "binary"}, the following model in considered for each SNP
  \deqn{ \mbox{logit}(P[Y=1| X, E, SNP, \omega])  = X\alpha + E\beta_E + G\beta_{SNP} + E*SNP\beta_{ExSNP} + \omega}{logit P(Y=1|X,E,SNP,omega)  = X alpha + E beta_E + SNP beta_SNP + E*SNP beta_ExSNP + omega}
  with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)}. 
  \eqn{E} the vector of the environmental variable, \eqn{SNP} is the genotype vector of the SNP, \eqn{K}{K} is a Genetic Relationship Matrix (GRM),
  \eqn{X} the covariable matrix. A score test is performed, independantly of the value of \code{test}.
  1 df score test is not availaible because it is very computationnaly heavy.
  All parameters under null model are estimated with the same procedure as in \code{\link{logistic.mm.aireml}}.
  It is possible to give a list in argument \code{K} for inclusion of several random effects in the model.
  If \code{p} is positive, the paramater \code{eigenK} needs to be given as well.
}

\value{
  A data frame, giving for each considered SNP, some of the following columns depending on the 
  values of \code{method} and \code{test}:
  \item{h2}{Estimated value of \eqn{\tau \over {\tau + \sigma^2}}{tau/(tau + sigma^2)}}
  \item{beta_E}{Estimated value of \eqn{\beta_E}{beta_E}}
  \item{sd_E}{Estimated standard deviation of the \eqn{\beta_E}{beta_E} estimation}
  \item{beta_SNP}{Estimated value of \eqn{\beta_{SNP}}{beta_SNP}}
  \item{sd_SNP}{Estimated standard deviation of the \eqn{\beta_{SNP}}{beta_SNP} estimation}
  \item{beta_ExSNP}{Estimated value of \eqn{\beta_{ExSNP}}{beta_ExNP}}
  \item{sd_ExSNP}{Estimated standard deviation of the \eqn{\beta_{ExSNP}}{beta_ExSNP} estimation}
  \item{cov_E_SNP}{Estimated covariance between the \eqn{\beta_E}{beta_E} and \eqn{\beta_{SNP}}{beta_SNP} estimations}
  \item{cov_E_ExSNP}{Estimated covariance between the \eqn{\beta_E}{beta_E} and \eqn{\beta_{SNP}}{beta_SNP} estimations}
  \item{cov_SNP_ExSNP}{Estimated covariance between the \eqn{\beta_{SNP}}{beta_SNP} and \eqn{\beta_{ExSNP}}{beta_ExSNP} estimations}  
  \item{Wald}{Wald statistic for each SNP}
  \item{score}{Score statistic for each SNP}
  \item{LRT}{Value of the Likelihood Ratio Test}
  \item{p}{The corresponding p-value}
}

\seealso{ \code{\link[gaston:qqplot.pvalues]{qqplot.pvalues}}, \code{\link[gaston:manhattan]{manhattan}},
          \code{\link[gaston:lmm.diago]{lmm.diago}}, \code{\link[gaston:lmm.aireml]{lmm.aireml}},
		  \code{\link[gaston:logistic.mm.aireml]{logistic.mm.aireml}}, \code{\link[stats:association.test]{association.test}} }

\examples{
\donttest{
# Load data of gaston package
data(TTN)
x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
standardize(x) <- "p"

# Simulate environmental variable
E <- rnorm(nrow(x), mean=5, sd=1)

# simulate quantitative phenotype with effect of SNP #631 and a polygenic component 
set.seed(1)
y <- E * (x \%*\% c(rep(0,630),0.5,rep(0,ncol(x)-631))) + rnorm(nrow(x))

# association test with linear model 
test <- genexE.association.test(x, y, E=E, method="lm", response = "quanti", df=1)

# a p-values qq plot
qqplot.pvalues(test$p)

# a small Manhattan plot 
# hihlighting the link between p-values and LD with SNP #631
lds <- LD(x, 631, c(1,ncol(x)))
manhattan(test, col = rgb(lds,0,0), pch = 20)

# use y to simulate a binary phenotype
y1 <- as.numeric(y > mean(y))

# logistic regression
t_binary <- genexE.association.test(x, y1, E=E, method = "lm", response = "binary", df=2)
# another small Manhattan plot for 2 df Wald test
t_binary$p <- t_binary$p
manhattan(t_binary, col = rgb(lds,0,0), pch = 20)
}}

\keyword{ Gene x Environment Iteracttion Association Test }
