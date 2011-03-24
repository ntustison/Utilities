###
### Function Code
###

'Bland.Altman' <- function(x,y,alpha=.05,rep.meas=FALSE,subject,...){
#**********************************************************************
#* Construct a Bland Altman Plot
#* 1. Set a few constants
#* 2. Calculate mean difference
#* 3. Calculate difference standard deviation
#* 4. Calculate upper and lower confidence limits
#* 5. Make Plot
#**********************************************************************

#*** 1. Set a few constants
  z <- qnorm(1-alpha/2)  ## value of z corresponding to alpha
  d <- x-y               ## pair-wise differences
  m <- (x+y)/2           ## pair-wise means

#*** 2. Calculate mean difference
  d.mn <- mean(d,na.rm=TRUE)

#*** 3. Calculate difference standard deviation
  if(rep.meas==FALSE){ d.sd=sqrt(var(d,na.rm=TRUE)) }
  else{

    #*** 3a. Ensure subject is a factor variable
    if(!is.factor(subject)) subject <- as.factor(subject)

    #*** 3b. Extract model information
    n <- length(levels(subject))      # Number of subjects
    model <- aov(d~subject)           # One way analysis of variance
    MSB <- anova(model)[[3]][1]       # Degrees of Freedom
    MSW <- anova(model)[[3]][2]       # Sums of Squares

    #*** 3c. Calculate number of complete pairs for each subject
    pairs <- NULL
    for(i in 1:length(levels(as.factor(subject)))){
      pairs[i] <- sum(is.na(d[subject==levels(subject)[i]])==FALSE)
    }
    Sig.dl <- (MSB-MSW)/((sum(pairs)^2-sum(pairs^2))/((n-1)*sum(pairs)))
    d.sd <- sqrt(Sig.dl+MSW)
  }

#*** 4. Calculate lower and upper confidence limits
  ucl <- d.mn+z*d.sd
  lcl <- d.mn-z*d.sd

#*** 5. Make Plot
  plot(m, d,abline(h=c(d.mn,ucl,lcl)),  ...)
  values <- round(cbind(lcl,d.mn,ucl),4)
  colnames(values) <- c("LCL","Mean","UCL")
  if(rep.meas==FALSE) Output <- list(limits=values,Var=d.sd^2)
    else Output <- list(limits=values,Var=Sig.dl)
  return(Output)
}


#
# Help File
#
#
#Bland Altman Plots
#
#Description:
#
#     Constructs a Bland-Altman Plot.
#
#Usage:
#
#      Bland.Altman(x,y,alpha=.05,rep.meas=FALSE,subject,...)
#
#Arguments:
#
#     x,y: vectors of values to be compared.
#
#   alpha: Significance level for determining confidence limits.
#          Defaults to 0.05
#
#rep.meas: Toggles if data provided should be considered as repeated
#          measures.  Defaults to 'FALSE'
#
# subject: Required if 'rep.meas=TRUE'.  A vector of the same length of
#          'x' and 'y' that  denotes which subject/group the measurement
#          belongs to.
#
#     ...: Other arguments to be passed to the 'plot' method.
#
#Details:
#
#     When 'rep.meas=TRUE', the confidence limits are calculated using a
#     method proposed by Bland and Altman. These limits are slightly
#     wider, allowing for the correlation within subjects/factors.  The
#     standard deviation used to compute these limits is:
#
#
#                sigma^2[d] = sigma^2[dI] + sigma^2[dw]
#
#
#     where  sigma^2[d]  is the variance of the differences, sigma^2[dI]
#      is the variance of the subjects and methods interaction, and
#     sigma^2[dw]  is the within subject variation. Estimates of these
#     values can be found with
#
#
#                            s^2[dw] = MSw
#
#
#
# s^2[dI] = (MSb - MSw) / ((sum(m[i])^2 - sum(m[i]^2)) /
#((n-1)*sum(m[i]))
#
#     )
#
#     Where MSb and MSw are the between and within subject variance of
#     the one way analysis of  variance and m[i] is the number of pairs
#     for the ith subject.  The sum of these two estimates provides the
#     estimate for  s^2[d]  .
#
#Value:
#
#  limits: A vector containing the Mean Bias and confidence limits.
#
#  Var.dl: The Variance of the Bias.  If 'rep.meas=TRUE', this is
#
#                               s^2[dI]
#
#          .
#
#Author(s):
#
#     Benjamin Nutter [hidden email]
#
#       Created:  December 2007
#
#References:
#
#     J Martin Bland and Douglass G Altman, "Measuring Agreement in
#     Method Comparison Studies", _Statistical Methods in Medical
#     Research_, 1999; 8: 135 - 160.
#
#     J Martin Bland and Doublas G. Altman, "Agreement Between Methods
#     of Measurement with Multiple Observations per Individual" _Journal
#     of Biopharmaceutical Statistics_ 17:571-582, 2007. (Corrects the
#     formula  given in the 1999 paper).
#
#     Burdick RK, Graybill FA. _Confidence Intervals on Variance
#     Components_. New York: Dekker, 1992.
#
#Examples:
#
#             observer1=rnorm(500,5,2)
#             observer2=rnorm(500,10,4)
#             ID=rep(1:50,10)
#
#             Bland.Altman(observer1,observer2)
#
#             Bland.Altman(observer1,observer2,rep.meas=TRUE,subject=ID)
