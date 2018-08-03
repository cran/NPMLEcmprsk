\name{NPMLEcmprsk}
\Rdversion{3.0}
\alias{NPMLEcmprsk}
\docType{package}
\title{Type-Specific Failure Rate and Hazard Rate on Competing Risks Data}
\description{
             Given a failure type, the function computes covariate-specific probability of failure 
             over time and covariate-specific conditional hazard rate based on possibly right-censored
             competing risk data. Specifically, it computes the non parametric maximum likelihood
             estimates of these quantities and their asymptotic variances in a semi parametric mixture
             model for competing risks data, as described in Chang et al. (2007a).
             }
\usage{
       NPMLEcmprsk(DATA,censoring.coding=0,alpha.stable.parameter=100, 
       beta.stable.parameter=100,initial.alpha=0,initial.beta=0, 
       threshold=0,iteration=5000)
       }
\arguments{
  \item{DATA}{
              An N-by-P matrix of data. Each of the N rows in the data matrix stores the data for one subject.
              Each of the P columns stores one variable. The first column stores the observation time, which is
              the minimum of time-to-failure and censoring time of a subject; the second column stores the observation
              time status, which indicates whether the observation time is a censoring time or not and the failure
              type of the subject when it is not a censoring time and is coded by an non-negative number for each
              observable time. The remaining columns store P-2 covariate values. The censoring status is coded zero.
              }
  \item{censoring.coding}{
                          The code represents the censoring status. Default value is zero.
                          }
  \item{alpha.stable.parameter,beta.stable.parameter}{
                                                     These parameters influence the convergence rate of the algorithm. More specifically, the parameter
                                                     determines the difference between the approximations of estimator in the algorithm; from the proposition
                                                     in Chang, et al. (2007a), the smaller parameter suggests the larger range between the approximations of
                                                     estimator. In our experiments, more covariates and more failure types usually have to provide larger stable
                                                     parameter for the convergence. Default value is 100.
                                                     }
  \item{initial.alpha,initial.beta}{
                                     The initial values of iterations in the algorithm. Default value is zero.
                                     }
  \item{threshold}{
                   The threshold for convergence in the algorithm. Default value is zero.
                   }
  \item{iteration}{
                   The number of iterations in the algorithm. Default value is 5,000.
                   }
}
\details{
         This function fits the semi-parametric mixture model for competing-risks data, in which proportional hazards models
         allowing covariates are specified for failure time models conditional on failure type and a multinomial model is specified
         for marginal distribution of failure type conditional on covariates, described in Chang et al. (2007a). This is a flexible
         model for competing risk data; in particular, it does not require the independence between covariates and censoring time.
         The algorithm is efficient for the computation of the non-parametric maximum likelihood estimates; it is also discussed in
         Chang et al. (2007b).
         }
\value{
       Returns a list with components
       \item{alpha}{
                    The coefficient on the probability of type-specific failures, defined by a logistic model, and it's a (P-1)-by-(K-1)
                    matrix, where K denote the number of failure types.
                    }
       \item{alpha.se,alpha.pvalue}{
                                    The standard error and p-value for alpha respectively.
                                    }              
       \item{alpha.95.lower.CI,alpha.95.upper.CI}{
                                                  The 95\% confidence interval for alpha.
                                                  }
       \item{beta}{
                   The relative risk coefficients in the hazard function for each type-specific time-to-failure and it's a (P-2)-by-K matrix.
                   }
       \item{beta.se,beta.pvalue}{
                                   The standard error and p-value for beta respectively.
                                   }
       \item{beta.95.lower.CI,beta.95.upper.CI}{
                                                The 95\% confidence interval for beta.
                                                }
       \item{Lambda}{
                     The cumulative baseline hazards.
                     }
       We will demonstrate a simulation study in the Examples section, in which included the plot of the probability of type-specific failures.
       }
\references{
            Chang, I. S. et al. (2007a) A fast algorithm for the nonparameteric maximum likelihood estimate in the Cox-gene model. Statist. Sinica 17, 841-856.\cr
            Chang, I. S. et al. (2007b) Non-parametric maximum-likelihood estimation in a semiparametric mixture model for competing-risks data. Scand. J. Statist. 34, 870-895. 
            }
\author{Chung-Hsing Chen <chchen@nhri.org.tw>}

\note{
      The missing value (e.g. NA) in the DATA is not allowed in this version.
      }

\seealso{
         none
         }

\examples{
\dontrun{
# setting the seed
set.seed(1)

# setting the sample size
N<-500

# setting the real parameters
alpha<-c(-2,5)
beta<-c(0.5,-0.5)
Lambda<-1/c(4,5)
Z<-rnorm(N,0,1)

# generating the data
W<-c((exp(cbind(1,Z) %*% alpha)/(1+exp(cbind(1,Z) %*% alpha))<runif(N))+1)
T<--log(matrix(runif(2*N),2,N))/Lambda*exp(-t(cbind(Z,Z))*beta)
survival.time<-sapply(1:N,function(i) T[W[i],i])

censoring.time<-runif(N,0.1,10)

temp<-survival.time<=censoring.time
X<-sapply(1:N,function(i) if(temp[i]) survival.time[i] else censoring.time[i])
delta<-sapply(1:N,function(i) if(temp[i]) W[i] else 0)

# estimating the parameters
data<-cbind(X,delta,Z)
result<-NPMLEcmprsk(data)
result

# plot probability of type-specific failures
logistic<-function(coef.alpha,covariates)
exp(cbind(1,covariates) %*% coef.alpha)/(1+exp(cbind(1,covariates) %*% coef.alpha))

range<-seq(min(Z),max(Z),0.1)

plot(range,logistic(result$coef.alpha,range),type="l"
,xlab="covariate",ylab="Case fatality rate")
lines(range,logistic(result$coef.alpha.95.lower.CI,range),lty=2)
lines(range,logistic(result$coef.alpha.95.upper.CI,range),lty=2)

points(Z,logistic(result$coef.alpha,Z),cex=0.7)
}
}


