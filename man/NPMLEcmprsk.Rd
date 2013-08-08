\name{NPMLEcmprsk}
\alias{NPMLEcmprsk}
\title{Non-parametric Maximum-Likelihood Estimation for Competing-Risks Data}
\description{
             The function computes non-parametric Maximum-Likelihood Estimate and its asymptotic variance 
             in a semiparametric mixture model for competing-risks data, as described in Chang et al. (2007b).
             }
\usage{NPMLEcmprsk(DATA, C, iteration)}
\arguments{
  \item{DATA}{
              The N-by-P matrix of data. There are N samples rows in matrix, with one sample in each 
              row. The P columns included the observable times which are time-to-event or censoring times, 
              the status which are coded by the non-negative number for each observable times, and the P-2 
              covariates. Note that the censoring status are coded by zero.
              }
  \item{C}{
           The parameter influence the convergence of the algorithm. The larger value suggests the faster 
           convergence. Default value is 50.
           }
           
  \item{iteration}{
                   The number of iteration for the algorithm. Default value is 300.
                   }
}
\details{
         Fits the semiparametric mixture model for competing-risks data, in which proportional hazards models 
         are specified for failure time models conditional on cause and a multinomial model is specified for 
         marginal distribution of cause conditional on covariates, described in Chang, I. S. et al. (2007b).
         The efficient algorithms for the computation of the non-parametric Maximum-Likelihood Estimation are 
         described in Change, I. S. et al. (2007a).
         }
\value{
       Returns a list with components
       \item{Lambda}{The cumulative hazard function. The format is a n-by-k matrix, where n is the number of 
                     event times and k is the type of status.}
       \item{alpha}{The coefficients on the probability of risk-sepcific failures, defined by a logistic model. 
                    The format is a (p-1)-by-(k-1) matrix, where p-1=(p-2)+1, the p-2 covariates and a constant.}
       \item{beta}{The coefficients on the survival function for time-to-event. The format is a (p-2)-by-k) matrix.}
       }
\references{
            Chang, I. S. et al. (2007a) A fast algorithm for the nonparameteric maximum likelihood estimate in the Cox-gene model. Statist. Sinica 17, 841-856.\cr
            Chang, I. S. et al. (2007b) Non-parametric maximum-likelihood estimation in a semiparametric mixture model for competing-risks data. Scand. J. Statist. 34, 870-895. 
            }
\author{Chung-Hsing Chen <chchen@nhri.org.tw>}

\note{
      The missing value (NA) in the DATA is not allowed in this version.
      }

\seealso{
         none
         }

\examples{
\dontrun{
# setting the seed
set.seed(1000)
                   
# setting the sample size
N=300
                   
# setting the real parameters
alpha=c(-2,5)
beta=c(0.5,-0.5)
Lambda=1/c(4,5)
Z=runif(N)
censoring.time=rexp(N,1/(50*(1-Z)))
                   
# generating the data
W=c((exp(cbind(1,Z) %*% alpha)/(1+exp(cbind(1,Z) %*% alpha))<runif(N))+1)
T=-log(matrix(runif(2*N),2,N))/Lambda*exp(-t(cbind(1,Z))*beta)
survival.time=sapply(1:N,function(i) T[W[i],i])
                   
temp=survival.time<=censoring.time
X=sapply(1:N,function(i) if(temp[i]) survival.time[i] else censoring.time[i])
delta=sapply(1:N,function(i) if(temp[i]) W[i] else 0)
                   
# estimating the parameters
data=cbind(X,delta,Z)
result=NP_competingRisk(data,0,300)
                   
result$alpha
result$beta
}
}

