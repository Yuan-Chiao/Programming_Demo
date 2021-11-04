###############################################################
#                                                             #
#                   Title: Injury Risk Functions              #
#                   Author: Yuan-Chiao Lu                     #
#                   Data sources: Simulated data              #
#                                                             #
###############################################################
rm(list=ls(all=TRUE))
GeneratedDataN <- 1

#############################################################################################
# Comparisons of 4 survival models: No Censor, Exact, Left-censored, and Interval-censored  #
#############################################################################################
Weibull.InterceptSE.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.InterceptCIRange.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.ScaleSE.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.ScaleCIRange.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.WeibullShapeSE.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.WeibullShapeCIRange.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.AlphaSE.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.AlphaCIRange.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime10SE.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime10CIRange.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime20SE.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime20CIRange.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime30SE.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime30CIRange.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime40SE.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime40CIRange.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.MedFailTimeSE.NoCensor.Successcount <- rep(0, GeneratedDataN)
Weibull.MedFailTimeCIRange.NoCensor.Successcount <- rep(0, GeneratedDataN)

Weibull.InterceptSE.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.InterceptCIRange.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.ScaleSE.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.ScaleCIRange.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.WeibullShapeSE.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.WeibullShapeCIRange.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.AlphaSE.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.AlphaCIRange.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime10SE.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime10CIRange.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime20SE.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime20CIRange.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime30SE.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime30CIRange.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime40SE.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime40CIRange.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.MedFailTimeSE.Exact.Successcount <- rep(0, GeneratedDataN)
Weibull.MedFailTimeCIRange.Exact.Successcount <- rep(0, GeneratedDataN)

Weibull.InterceptSE.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.InterceptCIRange.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.ScaleSE.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.ScaleCIRange.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.WeibullShapeSE.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.WeibullShapeCIRange.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.AlphaSE.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.AlphaCIRange.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime10SE.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime10CIRange.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime20SE.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime20CIRange.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime30SE.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime30CIRange.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime40SE.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime40CIRange.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.MedFailTimeSE.Left.Successcount <- rep(0, GeneratedDataN)
Weibull.MedFailTimeCIRange.Left.Successcount <- rep(0, GeneratedDataN)

Weibull.InterceptSE.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.InterceptCIRange.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.ScaleSE.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.ScaleCIRange.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.WeibullShapeSE.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.WeibullShapeCIRange.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.AlphaSE.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.AlphaCIRange.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime10SE.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime10CIRange.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime20SE.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime20CIRange.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime30SE.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime30CIRange.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime40SE.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.FailTime40CIRange.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.MedFailTimeSE.Interval.Successcount <- rep(0, GeneratedDataN)
Weibull.MedFailTimeCIRange.Interval.Successcount <- rep(0, GeneratedDataN)


for (OverallCounting in 1:GeneratedDataN){

################################################
#                                              #
#      Please Type in Initial Conditions       #
#                                              #
################################################
SampleN <- 20
Mu <- 350
Sigma <- 70
Min <- 250
Max <- 450
IntervalCensored_LowerBound <- 100
ProposedCTEStageN <- 1000
BootstrapN <- 2000
SurvCurvePointN <- 1000
################################################

###############################
#   Create Simulated data     #
###############################
RuptureSimulate <- Mu + Sigma*rnorm(1000)
AppliedSimulate <- runif(1000, Min, Max)
Rupture <- rep(0, SampleN)
Applied <- rep(0, SampleN)
Injury <- rep(0, SampleN)

k <- 1
for(i in 1:1000) {
	if(RuptureSimulate[i]<=AppliedSimulate[i]){
		Rupture[k] <- RuptureSimulate[i]
		Applied[k] <- AppliedSimulate[i]
		Injury[k] <- 1

		if (Rupture[k]<0){
			k <- k-1
		}
		k <- k+1
		if (k>(SampleN/2)){
			break 
		}
	}
}

k <- (SampleN/2)+1
for(i in 1:1000) {
	if(RuptureSimulate[i]>AppliedSimulate[i]){
		Rupture[k] <- RuptureSimulate[i]
		Applied[k] <- AppliedSimulate[i]
		Injury[k] <- 0
		k <- k+1
		if (k>SampleN){
			break
		}
	}
}


if (length(table(Applied))==SampleN){
	print("Good! No duplicated Applied Data!")
} else {
	print("Duplicated Applied Data!!!")
}

#Sort Applied Forces
n <- 1 # Because of "No duplicated Applied Data", so each row of dataset has only one trial (n=1). This is for Exact Logistic Regression.
Losses <- 1-Injury 
Deaths <- Injury
Stage <- Deaths/(Losses+Deaths)

###############################
#   Create Censored Data Set  #
###############################
Doubly_Lower <- rep(0,SampleN)
Doubly_Upper <- rep(0,SampleN)
Exact_Lower <- rep(0,SampleN)
Exact_Upper <- rep(0,SampleN)
Interval_Lower <- rep(0,SampleN)
Interval_Upper <- rep(0,SampleN)
Exact_Right_Data <- rep(0,SampleN)

for (i in 1:SampleN){
	if (Deaths[i]==0){
		Doubly_Lower[i] <- Applied[i]
		Doubly_Upper[i] <- NA
		Exact_Lower[i] <- Applied[i]
		Exact_Upper[i] <- NA
		Interval_Lower[i] <- Applied[i]
		Interval_Upper[i] <- NA
		Exact_Right_Data[i] <- Applied[i]
	} else{
		Doubly_Lower[i] <- NA
		Doubly_Upper[i] <- Applied[i]
		Exact_Lower[i] <- Rupture[i]
		Exact_Upper[i] <- Rupture[i]
		Interval_Lower[i] <- IntervalCensored_LowerBound
		Interval_Upper[i] <- Applied[i]
		Exact_Right_Data[i] <- Rupture[i]
	}
}

TestData <- cbind(Rupture, Applied, Losses, Deaths, Stage, n,
Doubly_Lower, Doubly_Upper, Exact_Lower, Exact_Upper, Interval_Lower, Interval_Upper, Exact_Right_Data)
SortTestData <- TestData[order(Applied) , ]

TestDataFrame <- as.data.frame(SortTestData, row.names = NULL)  #This is for Bootstrap
Applied <- SortTestData[,2]    #This is for Kaplan-Meier and Output Data
Losses <- SortTestData[,3]     #This is for Consistent Threshold Estimate
Deaths <- SortTestData[,4]     #This is for Consistent Threshold Estimate
InjuryData <- Deaths           #This is for Certainty Method
Stage <- SortTestData[,5]      #This is for Consistent Threshold Estimate
Deathinfo <- Deaths            #This is for Kaplan-Meier

write.table(SortTestData, file = "C:/Users/HEMI/Documents/GitHub/Programming_Demo/R/SortTestData.txt", col.names=NA)
SortTestData <- read.table("C:/Users/HEMI/Documents/GitHub/Programming_Demo/R/SortTestData.txt", head=T)

#################################################################
#                                                               #
#             Survival Data Analysis (Weibull)                  #
#                                                               #
#################################################################
####################################
#                                  #
#          All Exact Data          #
#         No Censored Data         #
#                                  #
####################################
library(splines)
library(survival)
wbs.simdata <- survreg(Surv(Rupture,n)~ 1, data=SortTestData, dist='weibull')
# dist = one of "extreme", "logistic", "gaussian", "weibull", "exponential", "rayleigh", "loggaussian", "lognormal", "loglogistic", "t"

multiplier <- qnorm(0.975, mean=0, sd=1, lower.tail = TRUE, log.p = FALSE)   # This is 1.96
InterceptMLE <- wbs.simdata$coeff
InterceptSE <- summary(wbs.simdata)$table[1,2]

ScaleMLE <- wbs.simdata$scale
ScaleSE <- sqrt((summary(wbs.simdata)$table[2,2]^2)*(wbs.simdata$scale^2))    # from Delta Method
logScaleMLE <- summary(wbs.simdata)$table[2,1]
logScaleSE <- summary(wbs.simdata)$table[2,2]
ScaleCI.Lower <- exp(logScaleMLE-multiplier*logScaleSE)
ScaleCI.Upper <- exp(logScaleMLE+multiplier*logScaleSE)

WeibullShapeMLE <- 1/wbs.simdata$scale
WeibullShapeSE <- sqrt((summary(wbs.simdata)$table[2,2]^2)*((-1/(exp(summary(wbs.simdata)$table[2,1])))^2))
WeibullShapeCI.Lower <- 1/exp(logScaleMLE+multiplier*logScaleSE)
WeibullShapeCI.Upper <- 1/exp(logScaleMLE-multiplier*logScaleSE)

CovInterceptScale <- ( (InterceptMLE^2)*(logScaleSE^2)+(2*logScaleMLE*wbs.simdata$var[1,2])-((InterceptMLE/ScaleMLE)^2)*(ScaleSE^2) )/(2*logScaleMLE*InterceptMLE/ScaleMLE)

AlphaMLE <- exp(-wbs.simdata$coeff/wbs.simdata$scale)
AlphaSE <- sqrt( (1/(ScaleMLE^2))*(AlphaMLE^2)*(InterceptSE^2)   +   ((InterceptMLE^2) / (ScaleMLE^4))*(AlphaMLE^2)*(ScaleSE^2)   +   (-2*InterceptMLE/(ScaleMLE^3))*(AlphaMLE^2)*CovInterceptScale )
AlphaCI.Lower <- exp(-(InterceptMLE+multiplier*InterceptSE)/ScaleCI.Lower) 
AlphaCI.Upper <- exp(-(InterceptMLE-multiplier*InterceptSE)/ScaleCI.Upper)

MedFailTimeMLE <- rep(0,SurvCurvePointN)
MedFailTimeSE <- rep(0,SurvCurvePointN)
MedFailTimeCI.Lower <- rep(0,SurvCurvePointN)
MedFailTimeCI.Upper <- rep(0,SurvCurvePointN)
Yindex <- rep(0,SurvCurvePointN)
CovWeibullShapeAlpha <- exp(-InterceptMLE/ScaleMLE)*(CovInterceptScale/(ScaleMLE^3)-InterceptMLE*(ScaleSE^2)/(ScaleMLE^4))
for (j in 1:(SurvCurvePointN-1)){
	Yindex[j+1] <- j
	MedFailTimeMLE[j+1] <- (log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaMLE)^(1/WeibullShapeMLE)
	MedFailTimeSE[j+1] <-   sqrt(   ((MedFailTimeMLE[j+1]*log(log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaMLE)*(-WeibullShapeMLE^(-2)))^2)*(WeibullShapeSE^2)  +
	(((1/WeibullShapeMLE)*((log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaMLE)^((1/WeibullShapeMLE)-1))*(-log(SurvCurvePointN/(SurvCurvePointN-j))*(AlphaMLE^(-2))))^2)*(AlphaSE^2)  +
	2*(MedFailTimeMLE[j+1]*log(log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaMLE)*(-WeibullShapeMLE^(-2)))*((1/WeibullShapeMLE)*((log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaMLE)^((1/WeibullShapeMLE)-1))*(-log(SurvCurvePointN/(SurvCurvePointN-j))*(AlphaMLE^(-2))))*CovWeibullShapeAlpha   )
	#MedFailTimeCI.Lower[j+1] <- (log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaCI.Upper)^(1/WeibullShapeCI.Lower)
	#MedFailTimeCI.Upper[j+1] <- (log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaCI.Lower)^(1/WeibullShapeCI.Upper)
	MedFailTimeCI.Lower[j+1] <- MedFailTimeMLE[j+1]-multiplier*MedFailTimeSE[j+1]
	MedFailTimeCI.Upper[j+1] <- MedFailTimeMLE[j+1]+multiplier*MedFailTimeSE[j+1]
}
InjuryTimeData.NoCensor <- cbind(Yindex,MedFailTimeMLE,MedFailTimeSE,MedFailTimeCI.Lower,MedFailTimeCI.Upper)
write.csv(InjuryTimeData.NoCensor, file = "C:/Users/HEMI/Documents/GitHub/Programming_Demo/R/InjuryTimeData.NoCensor.csv")

InterceptOutput.NoCensor <- cbind(InterceptMLE, InterceptSE, InterceptMLE-multiplier*InterceptSE, InterceptMLE+multiplier*InterceptSE, 2*multiplier*InterceptSE, (multiplier*InterceptSE)/(multiplier*InterceptSE))

ScaleOutput.NoCensor <- cbind(ScaleMLE, ScaleSE, ScaleCI.Lower, ScaleCI.Upper, ScaleCI.Upper-ScaleCI.Lower, (ScaleCI.Upper-ScaleMLE)/(ScaleMLE-ScaleCI.Lower))

WeibullShapeCIRange <- WeibullShapeCI.Upper-WeibullShapeCI.Lower
WeibullShapeCIAsymmetryIndex <- (WeibullShapeCI.Upper-WeibullShapeMLE)/(WeibullShapeMLE-WeibullShapeCI.Lower)
WeibullShapeOutput.NoCensor <- cbind(WeibullShapeMLE, WeibullShapeSE, WeibullShapeCI.Lower, WeibullShapeCI.Upper, WeibullShapeCIRange, WeibullShapeCIAsymmetryIndex)

AlphaCIRange <- AlphaCI.Upper-AlphaCI.Lower
AlphaCIAsymmetryIndex <- (AlphaCI.Upper-AlphaMLE)/(AlphaMLE-AlphaCI.Lower)
AlphaOutput.NoCensor <- cbind(AlphaMLE, AlphaSE, AlphaCI.Lower, AlphaCI.Upper, AlphaCIRange, AlphaCIAsymmetryIndex)

FailTime10CIRange <- MedFailTimeCI.Upper[SurvCurvePointN/10+1]-MedFailTimeCI.Lower[SurvCurvePointN/10+1]
FailTime10CIAsymmetryIndex <- (MedFailTimeCI.Upper[SurvCurvePointN/10+1]-MedFailTimeMLE[SurvCurvePointN/10+1])/(MedFailTimeMLE[SurvCurvePointN/10+1]-MedFailTimeCI.Lower[SurvCurvePointN/10+1])
FailTime10Output.NoCensor <- cbind(MedFailTimeMLE[SurvCurvePointN/10+1],MedFailTimeSE[SurvCurvePointN/10+1],MedFailTimeCI.Lower[SurvCurvePointN/10+1],MedFailTimeCI.Upper[SurvCurvePointN/10+1],FailTime10CIRange, FailTime10CIAsymmetryIndex)
FailTime20CIRange <- MedFailTimeCI.Upper[SurvCurvePointN/5+1]-MedFailTimeCI.Lower[SurvCurvePointN/5+1]
FailTime20CIAsymmetryIndex <- (MedFailTimeCI.Upper[SurvCurvePointN/5+1]-MedFailTimeMLE[SurvCurvePointN/5+1])/(MedFailTimeMLE[SurvCurvePointN/5+1]-MedFailTimeCI.Lower[SurvCurvePointN/5+1])
FailTime20Output.NoCensor <- cbind(MedFailTimeMLE[SurvCurvePointN/5+1],MedFailTimeSE[SurvCurvePointN/5+1],MedFailTimeCI.Lower[SurvCurvePointN/5+1],MedFailTimeCI.Upper[SurvCurvePointN/5+1],FailTime20CIRange, FailTime20CIAsymmetryIndex)
FailTime30CIRange <- MedFailTimeCI.Upper[SurvCurvePointN/(10/3)+1]-MedFailTimeCI.Lower[SurvCurvePointN/(10/3)+1]
FailTime30CIAsymmetryIndex <- (MedFailTimeCI.Upper[SurvCurvePointN/(10/3)+1]-MedFailTimeMLE[SurvCurvePointN/(10/3)+1])/(MedFailTimeMLE[SurvCurvePointN/(10/3)+1]-MedFailTimeCI.Lower[SurvCurvePointN/(10/3)+1])
FailTime30Output.NoCensor <- cbind(MedFailTimeMLE[SurvCurvePointN/(10/3)+1],MedFailTimeSE[SurvCurvePointN/(10/3)+1],MedFailTimeCI.Lower[SurvCurvePointN/(10/3)+1],MedFailTimeCI.Upper[SurvCurvePointN/(10/3)+1],FailTime30CIRange, FailTime30CIAsymmetryIndex)
FailTime40CIRange <- MedFailTimeCI.Upper[SurvCurvePointN/2.5+1]-MedFailTimeCI.Lower[SurvCurvePointN/2.5+1]
FailTime40CIAsymmetryIndex <- (MedFailTimeCI.Upper[SurvCurvePointN/2.5+1]-MedFailTimeMLE[SurvCurvePointN/2.5+1])/(MedFailTimeMLE[SurvCurvePointN/2.5+1]-MedFailTimeCI.Lower[SurvCurvePointN/2.5+1])
FailTime40Output.NoCensor <- cbind(MedFailTimeMLE[SurvCurvePointN/2.5+1],MedFailTimeSE[SurvCurvePointN/2.5+1],MedFailTimeCI.Lower[SurvCurvePointN/2.5+1],MedFailTimeCI.Upper[SurvCurvePointN/2.5+1],FailTime40CIRange, FailTime40CIAsymmetryIndex)
MedFailTimeCIRange <- MedFailTimeCI.Upper[SurvCurvePointN/2+1]-MedFailTimeCI.Lower[SurvCurvePointN/2+1]
MedFailTimeCIAsymmetryIndex <- (MedFailTimeCI.Upper[SurvCurvePointN/2+1]-MedFailTimeMLE[SurvCurvePointN/2+1])/(MedFailTimeMLE[SurvCurvePointN/2+1]-MedFailTimeCI.Lower[SurvCurvePointN/2+1])
MedFailTimeOutput.NoCensor <- cbind(MedFailTimeMLE[SurvCurvePointN/2+1],MedFailTimeSE[SurvCurvePointN/2+1],MedFailTimeCI.Lower[SurvCurvePointN/2+1],MedFailTimeCI.Upper[SurvCurvePointN/2+1],MedFailTimeCIRange,MedFailTimeCIAsymmetryIndex)

MeanInjuryProb.Weibull.NoCensor <- rep(0,Mu*2+1)
SEInjuryProb.Weibull.NoCensor <- rep(0,Mu*2+1)
LowerInjuryProb.Weibull.NoCensor <- rep(0,Mu*2+1)
UpperInjuryProb.Weibull.NoCensor <- rep(0,Mu*2+1)
Xindex <- rep(0,Mu*2+1)
CovWeibullShapeAlpha <- exp(-InterceptMLE/ScaleMLE)*(CovInterceptScale/(ScaleMLE^3)-InterceptMLE*(ScaleSE^2)/(ScaleMLE^4))
CovMatrix <- rbind(c(WeibullShapeSE^2,CovWeibullShapeAlpha),c(CovWeibullShapeAlpha,AlphaSE^2))
for (j in 1:(Mu*2)){
	Xindex[j+1] <- j
	MeanInjuryProb.Weibull.NoCensor[j+1] <- 1-exp(-AlphaMLE*(j^WeibullShapeMLE))
	DiffWeibullShapeAlpha <- c( AlphaMLE*(j^WeibullShapeMLE)*log(j)*exp(-AlphaMLE*(j^WeibullShapeMLE)), (j^WeibullShapeMLE)*exp(-AlphaMLE*(j^WeibullShapeMLE)) )
	SEInjuryProb.Weibull.NoCensor[j+1] <- sqrt( t(DiffWeibullShapeAlpha) %*% CovMatrix %*% DiffWeibullShapeAlpha )
	LowerInjuryProb.Weibull.NoCensor[j+1] <- MeanInjuryProb.Weibull.NoCensor[j+1]-1.96*SEInjuryProb.Weibull.NoCensor[j+1]
	UpperInjuryProb.Weibull.NoCensor[j+1] <- MeanInjuryProb.Weibull.NoCensor[j+1]+1.96*SEInjuryProb.Weibull.NoCensor[j+1]
}
InjuryProbData <- cbind(Xindex,MeanInjuryProb.Weibull.NoCensor,LowerInjuryProb.Weibull.NoCensor,UpperInjuryProb.Weibull.NoCensor)
write.csv(InjuryProbData, file = "C:/Users/HEMI/Documents/GitHub/Programming_Demo/R/InjuryProbData.csv")


####################################
#                                  #
#        Injury: Exact              #
#    Non-Injury: Right Censored    #
#                                  #
####################################
wbs.Exact_Right <- survreg(Surv(Exact_Right_Data,Deaths)~ 1, data=SortTestData, dist='weibull')
# Note: All other codes are similar to "All Exact Data" model. Skip for demonstration purpose.

####################################
#                                  #
#        Injury: Left Censored     #
#    Non-Injury: Right Censored    #
#                                  #
####################################
wbs.LeftCensored <- survreg(Surv(Doubly_Lower,Doubly_Upper,type="interval2")~ 1, data=SortTestData, dist='weibull')
# Note: All other codes are similar to "All Exact Data" model. Skip for demonstration purpose.

####################################
#                                  #
#        Injury: Interval Censored #
#    Non-Injury: Right Censored    #
#                                  #
####################################
wbs.interval <- survreg(Surv(Interval_Lower,Interval_Upper,type="interval2")~ 1, data=SortTestData, dist='weibull')
# Note: All other codes are similar to "All Exact Data" model. Skip for demonstration purpose.

#################################################################
#                                                               #
#    Bootstrap 95% CI of parameters Alpha and Weibull Shape     #
#                                                               #
#################################################################
RuptureBoot <- rep(0, SampleN)
DeathsBoot <- rep(0, SampleN)
n <- rep(1, SampleN)
Exact_Right_DataBoot <- rep(0, SampleN)
Doubly_Lower_InfinityBoot <- rep(0, SampleN)
Doubly_Upper_InfinityBoot <- rep(0, SampleN)
Interval_LowerBoot <- rep(0, SampleN)
Interval_Upper_InfinityBoot <- rep(0, SampleN)

Weibull.NoCensor.coeff.boot <- rep(0, BootstrapN)
Weibull.NoCensor.scale.boot <- rep(0, BootstrapN)
Weibull.Exact.coeff.boot <- rep(0, BootstrapN)
Weibull.Exact.scale.boot <- rep(0, BootstrapN)
Weibull.LeftCensored.coeff.boot <- rep(0, BootstrapN)
Weibull.LeftCensored.scale.boot <- rep(0, BootstrapN)
Weibull.IntervalCensored.coeff.boot <- rep(0, BootstrapN)
Weibull.IntervalCensored.scale.boot <- rep(0, BootstrapN)

k <- 1
Testindex1 <- 0
Testindex2 <- 0
Testindex3 <- 0
Testindex4 <- 0
Testindex5 <- 0
Testindex6 <- 0
for (j in 1:(BootstrapN*3)){    #BootstrapN*3 is a randomly assigned number, which is greater than BootstrapN
	BootIndex <- 1
	SkipIndex <- 0
	for (i in 1:SampleN) {
		RandomN <- sample(1:SampleN, 1)
     		RuptureBoot[i] <- TestDataFrame$Rupture[RandomN]

		RandomN <- sample(1:(SampleN/2),1)
		if (i <= (SampleN/2)){
			DeathsBoot[i] <- 1                                             #resample the injury data
			Exact_Right_DataBoot[i] <- TestData[RandomN,1]
			Doubly_Lower_InfinityBoot[i] <- NA
			Doubly_Upper_InfinityBoot[i] <- TestData[RandomN,2]
			Interval_LowerBoot[i] <- IntervalCensored_LowerBound
			Interval_Upper_InfinityBoot[i] <- TestData[RandomN,2]
			} else {
			DeathsBoot[i] <- 0                                             #resample the non-injury data
			Exact_Right_DataBoot[i] <- TestData[RandomN+(SampleN/2),2]
			Doubly_Lower_InfinityBoot[i] <- TestData[RandomN+(SampleN/2),2]
			Doubly_Upper_InfinityBoot[i] <- NA
			Interval_LowerBoot[i] <- TestData[RandomN+(SampleN/2),2]
			Interval_Upper_InfinityBoot[i] <- NA
		}
	}
	TestDataboot <- cbind(RuptureBoot, n, DeathsBoot,Exact_Right_DataBoot,Doubly_Lower_InfinityBoot,
	Doubly_Upper_InfinityBoot,Interval_LowerBoot,Interval_Upper_InfinityBoot)

	TestDataFrameboot <- as.data.frame(TestDataboot, row.names = NULL)  #This is for Bootstrap

	if (BootIndex==1){

		Weibull.NoCensor.boot <- survreg(Surv(RuptureBoot,n)~ 1, data=TestDataFrameboot , dist='weibull')
		Weibull.NoCensor.coeff.boot[k] <- Weibull.NoCensor.boot$coeff
		Weibull.NoCensor.scale.boot[k] <- Weibull.NoCensor.boot$scale

		Weibull.Exact.boot <- survreg(Surv(Exact_Right_DataBoot,DeathsBoot)~ 1, data=TestDataFrameboot , dist='weibull')
		Weibull.Exact.coeff.boot[k] <- Weibull.Exact.boot$coeff
		Weibull.Exact.scale.boot[k] <- Weibull.Exact.boot$scale

		Weibull.LeftCensored.boot <- survreg(Surv(Doubly_Lower_InfinityBoot,Doubly_Upper_InfinityBoot,type="interval2")~ 1, data=TestDataFrameboot , dist='weibull')
		Weibull.LeftCensored.coeff.boot[k] <- Weibull.LeftCensored.boot$coeff
		Weibull.LeftCensored.scale.boot[k] <- Weibull.LeftCensored.boot$scale

		Weibull.IntervalCensored.boot <- survreg(Surv(Interval_LowerBoot,Interval_Upper_InfinityBoot,type="interval2")~ 1, data=TestDataFrameboot , dist='weibull')
		Weibull.IntervalCensored.coeff.boot[k] <- Weibull.IntervalCensored.boot$coeff
		Weibull.IntervalCensored.scale.boot[k] <- Weibull.IntervalCensored.boot$scale


		if (Weibull.LeftCensored.boot$iter >= 30){         #If iterations=30, it shows "Ran out of iterations and did not converge".
			SkipIndex <- 1
			Testindex1 <- Testindex1+1
		}
		if (Weibull.LeftCensored.boot$scale=="NaN") {
			SkipIndex <- 1
			Testindex3 <- Testindex3+1	
		} else {
			if (Weibull.LeftCensored.boot$scale >= 20){
				SkipIndex <- 1
			}
			if ((exp(-Weibull.LeftCensored.boot$coeff/Weibull.LeftCensored.boot$scale)) < 1e-300){
				SkipIndex <- 1
				Testindex2 <- Testindex2+1
			}
		}
		if (Weibull.IntervalCensored.boot$iter >= 30){
			SkipIndex <- 1
			Testindex4 <- Testindex4+1
		}
		if (Weibull.IntervalCensored.boot$scale=="NaN") {
			SkipIndex <- 1
			Testindex6 <- Testindex6+1		
		} else {
			if (Weibull.IntervalCensored.boot$scale >= 20){
				SkipIndex <- 1
			}
			if ((exp(-Weibull.IntervalCensored.boot$coeff/Weibull.IntervalCensored.boot$scale)) < 1e-300){
				SkipIndex <- 1
				Testindex5 <- Testindex5+1
			}
		}
		if (SkipIndex==1){
			k <- k-1
		}

	k <- k+1
	}
	if (k==(BootstrapN+1)){
		break
	}
}

Testindex1
Testindex2
Testindex3
Testindex4
Testindex5
Testindex6
k


###############################
#         Injury: Exact       #
###############################
Weibull.Exact.Alpha.boot <- exp(-Weibull.Exact.coeff.boot/Weibull.Exact.scale.boot)
Weibull.Exact.Alpha.boot.SurvCurve <- Weibull.Exact.Alpha.boot
Weibull.Exact.WeibullShape.boot <- 1/Weibull.Exact.scale.boot
Weibull.Exact.WeibullShape.boot.SurvCurve <- Weibull.Exact.WeibullShape.boot


InterceptMLEboot <- mean(Weibull.Exact.coeff.boot)
InterceptSEboot <- sd(Weibull.Exact.coeff.boot)
Weibull.Exact.coeff.boot <- sort(Weibull.Exact.coeff.boot)
#95% CI for Intercept is (Lower Bound): Weibull.Exact.coeff.boot[BootstrapN*0.05/2]
#95% CI for Intercept is (Upper Bound): Weibull.Exact.coeff.boot[BootstrapN*(1-0.05/2)]

ScaleMLEboot <- mean(Weibull.Exact.scale.boot)
ScaleSEboot <- sd(Weibull.Exact.scale.boot)
Weibull.Exact.scale.boot <- sort(Weibull.Exact.scale.boot)
#95% CI for Scale is (Lower Bound): Weibull.Exact.scale.boot[BootstrapN*0.05/2]
#95% CI for Scale is (Upper Bound): Weibull.Exact.scale.boot[BootstrapN*(1-0.05/2)]

WeibullShapeMLEboot <- 1/ScaleMLEboot
WeibullShapeSEboot <- sqrt((ScaleSEboot^2)*((-1/(ScaleMLEboot^2))^2))
Weibull.Exact.WeibullShape.boot <- sort(Weibull.Exact.WeibullShape.boot)
WeibullShapeCIboot.Lower <- Weibull.Exact.WeibullShape.boot[BootstrapN*0.05/2]
WeibullShapeCIboot.Upper <- Weibull.Exact.WeibullShape.boot[BootstrapN*(1-0.05/2)]

AlphaMLEboot <- exp(-InterceptMLEboot/ScaleMLEboot)
AlphaSEboot <- sqrt( (1/(ScaleMLEboot^2))*(AlphaMLEboot^2)*(InterceptSEboot^2)   +   ((InterceptMLEboot^2) / (ScaleMLEboot^4))*(AlphaMLEboot^2)*(ScaleSEboot^2)   +   (-2*InterceptMLEboot/(ScaleMLEboot^3))*(AlphaMLEboot^2)*cov(Weibull.Exact.coeff.boot,Weibull.Exact.scale.boot) )
Weibull.Exact.Alpha.boot <- sort(Weibull.Exact.Alpha.boot)
AlphaCIboot.Lower <- Weibull.Exact.Alpha.boot[BootstrapN*0.05/2]
AlphaCIboot.Upper <- Weibull.Exact.Alpha.boot[BootstrapN*(1-0.05/2)]



MedFailTimeMLEboot <- rep(0,SurvCurvePointN)
MedFailTimeSEboot <- rep(0,SurvCurvePointN)
MedFailTimeCIboot.Lower <- rep(0,SurvCurvePointN)
MedFailTimeCIboot.Upper <- rep(0,SurvCurvePointN)
Yindex <- rep(0,SurvCurvePointN)
CovWeibullShapeAlphaboot <- exp(-InterceptMLEboot/ScaleMLEboot)*(cov(Weibull.Exact.coeff.boot,Weibull.Exact.scale.boot)/(ScaleMLEboot^3)-InterceptMLEboot*(ScaleSEboot^2)/(ScaleMLEboot^4))
for (j in 1:(SurvCurvePointN-1)){
	Yindex[j+1] <- j
	#MedFailTimeMLEboot[j+1] <- (log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaMLEboot)^(1/WeibullShapeMLEboot)
	#MedFailTimeSEboot[j+1] <- sqrt(   ((MedFailTimeMLEboot[j+1]*log(log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaMLEboot)*(-WeibullShapeMLEboot^(-2)))^2)*(WeibullShapeSEboot^2)  +
	#(((1/WeibullShapeMLEboot)*((log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaMLEboot)^((1/WeibullShapeMLEboot)-1))*(-log(SurvCurvePointN/(SurvCurvePointN-j))*(AlphaMLEboot^(-2))))^2)*(AlphaSEboot^2)  +
	#2*(MedFailTimeMLEboot[j+1]*log(log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaMLEboot)*(-WeibullShapeMLEboot^(-2)))*((1/WeibullShapeMLEboot)*((log(SurvCurvePointN/(SurvCurvePointN-j))/AlphaMLEboot)^((1/WeibullShapeMLEboot)-1))*(-log(SurvCurvePointN/(SurvCurvePointN-j))*(AlphaMLEboot^(-2))))*CovWeibullShapeAlphaboot   )
	Weibull.Exact.MedFailTime.boot <- (log(SurvCurvePointN/(SurvCurvePointN-j))/Weibull.Exact.Alpha.boot.SurvCurve)^(1/Weibull.Exact.WeibullShape.boot.SurvCurve)
	MedFailTimeMLEboot[j+1] <- mean(Weibull.Exact.MedFailTime.boot)
	MedFailTimeSEboot[j+1] <- sd(Weibull.Exact.MedFailTime.boot)
	Weibull.Exact.MedFailTime.boot <- sort(Weibull.Exact.MedFailTime.boot)
	MedFailTimeCIboot.Lower[j+1] <- Weibull.Exact.MedFailTime.boot[BootstrapN*0.05/2]
	MedFailTimeCIboot.Upper[j+1] <- Weibull.Exact.MedFailTime.boot[BootstrapN*(1-0.05/2)]
}
InjuryTimeDataboot.Exact <- cbind(Yindex,MedFailTimeMLEboot,MedFailTimeSEboot,MedFailTimeCIboot.Lower,MedFailTimeCIboot.Upper)
write.csv(InjuryTimeDataboot.Exact, file = "C:/Users/HEMI/Documents/GitHub/Programming_Demo/R/InjuryTimeDataboot.Exact.csv")

WeibullShapeCIbootRange <- WeibullShapeCIboot.Upper-WeibullShapeCIboot.Lower
WeibullShapeCIbootAsymmetryIndex <- (WeibullShapeCIboot.Upper-WeibullShapeMLEboot)/(WeibullShapeMLEboot-WeibullShapeCIboot.Lower)
WeibullShapebootOutput.Exact <- cbind(WeibullShapeMLEboot , WeibullShapeSEboot, WeibullShapeCIboot.Lower, WeibullShapeCIboot.Upper, WeibullShapeCIbootRange, WeibullShapeCIbootAsymmetryIndex)

AlphaCIbootRange <- AlphaCIboot.Upper-AlphaCIboot.Lower
AlphaCIbootAsymmetryIndex <- (AlphaCIboot.Upper-AlphaMLEboot)/(AlphaMLEboot-AlphaCIboot.Lower)
AlphabootOutput.Exact <- cbind(AlphaMLEboot, AlphaSEboot, AlphaCIboot.Lower, AlphaCIboot.Upper, AlphaCIbootRange, AlphaCIbootAsymmetryIndex)

FailTime10CIbootRange <- MedFailTimeCIboot.Upper[SurvCurvePointN/10+1]-MedFailTimeCIboot.Lower[SurvCurvePointN/10+1]
FailTime10CIbootAsymmetryIndex <- (MedFailTimeCIboot.Upper[SurvCurvePointN/10+1]-MedFailTimeMLEboot[SurvCurvePointN/10+1])/(MedFailTimeMLEboot[SurvCurvePointN/10+1]-MedFailTimeCIboot.Lower[SurvCurvePointN/10+1])
FailTime10bootOutput.Exact <- cbind(MedFailTimeMLEboot[SurvCurvePointN/10+1],MedFailTimeSEboot[SurvCurvePointN/10+1],MedFailTimeCIboot.Lower[SurvCurvePointN/10+1],MedFailTimeCIboot.Upper[SurvCurvePointN/10+1], FailTime10CIbootRange,FailTime10CIbootAsymmetryIndex)
FailTime20CIbootRange <- MedFailTimeCIboot.Upper[SurvCurvePointN/5+1]-MedFailTimeCIboot.Lower[SurvCurvePointN/5+1]
FailTime20CIbootAsymmetryIndex <- (MedFailTimeCIboot.Upper[SurvCurvePointN/5+1]-MedFailTimeMLEboot[SurvCurvePointN/5+1])/(MedFailTimeMLEboot[SurvCurvePointN/5+1]-MedFailTimeCIboot.Lower[SurvCurvePointN/5+1])
FailTime20bootOutput.Exact <- cbind(MedFailTimeMLEboot[SurvCurvePointN/5+1],MedFailTimeSEboot[SurvCurvePointN/5+1],MedFailTimeCIboot.Lower[SurvCurvePointN/5+1],MedFailTimeCIboot.Upper[SurvCurvePointN/5+1], FailTime20CIbootRange,FailTime20CIbootAsymmetryIndex)
FailTime30CIbootRange <- MedFailTimeCIboot.Upper[SurvCurvePointN/(10/3)+1]-MedFailTimeCIboot.Lower[SurvCurvePointN/(10/3)+1]
FailTime30CIbootAsymmetryIndex <- (MedFailTimeCIboot.Upper[SurvCurvePointN/(10/3)+1]-MedFailTimeMLEboot[SurvCurvePointN/(10/3)+1])/(MedFailTimeMLEboot[SurvCurvePointN/(10/3)+1]-MedFailTimeCIboot.Lower[SurvCurvePointN/(10/3)+1])
FailTime30bootOutput.Exact <- cbind(MedFailTimeMLEboot[SurvCurvePointN/(10/3)+1],MedFailTimeSEboot[SurvCurvePointN/(10/3)+1],MedFailTimeCIboot.Lower[SurvCurvePointN/(10/3)+1],MedFailTimeCIboot.Upper[SurvCurvePointN/(10/3)+1], FailTime30CIbootRange,FailTime30CIbootAsymmetryIndex)
FailTime40CIbootRange <- MedFailTimeCIboot.Upper[SurvCurvePointN/2.5+1]-MedFailTimeCIboot.Lower[SurvCurvePointN/2.5+1]
FailTime40CIbootAsymmetryIndex <- (MedFailTimeCIboot.Upper[SurvCurvePointN/2.5+1]-MedFailTimeMLEboot[SurvCurvePointN/2.5+1])/(MedFailTimeMLEboot[SurvCurvePointN/2.5+1]-MedFailTimeCIboot.Lower[SurvCurvePointN/2.5+1])
FailTime40bootOutput.Exact <- cbind(MedFailTimeMLEboot[SurvCurvePointN/2.5+1],MedFailTimeSEboot[SurvCurvePointN/2.5+1],MedFailTimeCIboot.Lower[SurvCurvePointN/2.5+1],MedFailTimeCIboot.Upper[SurvCurvePointN/2.5+1], FailTime40CIbootRange,FailTime40CIbootAsymmetryIndex)
MedFailTimeCIbootRange <- MedFailTimeCIboot.Upper[SurvCurvePointN/2+1]-MedFailTimeCIboot.Lower[SurvCurvePointN/2+1]
MedFailTimeCIbootAsymmetryIndex <- (MedFailTimeCIboot.Upper[SurvCurvePointN/2+1]-MedFailTimeMLEboot[SurvCurvePointN/2+1])/(MedFailTimeMLEboot[SurvCurvePointN/2+1]-MedFailTimeCIboot.Lower[SurvCurvePointN/2+1])
MedFailTimebootOutput.Exact <- cbind(MedFailTimeMLEboot[SurvCurvePointN/2+1],MedFailTimeSEboot[SurvCurvePointN/2+1],MedFailTimeCIboot.Lower[SurvCurvePointN/2+1],MedFailTimeCIboot.Upper[SurvCurvePointN/2+1],MedFailTimeCIbootRange,MedFailTimeCIbootAsymmetryIndex)



if (InterceptSEboot < InterceptOutput.Exact[2]){
Weibull.InterceptSE.Exact.Successcount[OverallCounting] <- 1
}
if ((Weibull.Exact.coeff.boot[BootstrapN*(1-0.05/2)]-Weibull.Exact.coeff.boot[BootstrapN*0.05/2]) < InterceptOutput.Exact[5]){
Weibull.InterceptCIRange.Exact.Successcount[OverallCounting] <- 1
}
if (ScaleSEboot < ScaleOutput.Exact[2]){
Weibull.ScaleSE.Exact.Successcount[OverallCounting] <- 1
}
if ((Weibull.Exact.scale.boot[BootstrapN*(1-0.05/2)]-Weibull.Exact.scale.boot[BootstrapN*0.05/2]) < ScaleOutput.Exact[5]){
Weibull.ScaleCIRange.Exact.Successcount[OverallCounting] <- 1
}

if (WeibullShapeSEboot < WeibullShapeOutput.Exact[2]){
Weibull.WeibullShapeSE.Exact.Successcount[OverallCounting] <- 1
}
if (WeibullShapeCIbootRange < WeibullShapeOutput.Exact[5]){
Weibull.WeibullShapeCIRange.Exact.Successcount[OverallCounting] <- 1
}
if (AlphaSEboot < AlphaOutput.Exact[2]){
Weibull.AlphaSE.Exact.Successcount[OverallCounting] <- 1
}
if (AlphaCIbootRange < AlphaOutput.Exact[5]){
Weibull.AlphaCIRange.Exact.Successcount[OverallCounting] <- 1
}


if (FailTime10bootOutput.Exact[2] < FailTime10Output.Exact[2]){
Weibull.FailTime10SE.Exact.Successcount[OverallCounting] <- 1
}
if (FailTime10bootOutput.Exact[5] < FailTime10Output.Exact[5]){
Weibull.FailTime10CIRange.Exact.Successcount[OverallCounting] <- 1
}
if (FailTime20bootOutput.Exact[2] < FailTime20Output.Exact[2]){
Weibull.FailTime20SE.Exact.Successcount[OverallCounting] <- 1
}
if (FailTime20bootOutput.Exact[5] < FailTime20Output.Exact[5]){
Weibull.FailTime20CIRange.Exact.Successcount[OverallCounting] <- 1
}
if (FailTime30bootOutput.Exact[2] < FailTime30Output.Exact[2]){
Weibull.FailTime30SE.Exact.Successcount[OverallCounting] <- 1
}
if (FailTime30bootOutput.Exact[5] < FailTime30Output.Exact[5]){
Weibull.FailTime30CIRange.Exact.Successcount[OverallCounting] <- 1
}
if (FailTime40bootOutput.Exact[2] < FailTime40Output.Exact[2]){
Weibull.FailTime40SE.Exact.Successcount[OverallCounting] <- 1
}
if (FailTime40bootOutput.Exact[5] < FailTime40Output.Exact[5]){
Weibull.FailTime40CIRange.Exact.Successcount[OverallCounting] <- 1
}
if (MedFailTimebootOutput.Exact[2] < MedFailTimeOutput.Exact[2]){
Weibull.MedFailTimeSE.Exact.Successcount[OverallCounting] <- 1
}
if (MedFailTimebootOutput.Exact[5] < MedFailTimeOutput.Exact[5]){
Weibull.MedFailTimeCIRange.Exact.Successcount[OverallCounting] <- 1
}

####################################
#                                  #
#        Injury: Exact              #
#    Non-Injury: Right Censored    #
#                                  #
####################################
# Note: Codes are similar to "All Exact Data" model. Skip for demonstration purpose.

####################################
#                                  #
#        Injury: Left Censored     #
#    Non-Injury: Right Censored    #
#                                  #
####################################
# Note: Codes are similar to "All Exact Data" model. Skip for demonstration purpose.

####################################
#                                  #
#        Injury: Interval Censored #
#    Non-Injury: Right Censored    #
#                                  #
####################################
# Note: Codes are similar to "All Exact Data" model. Skip for demonstration purpose.

###############################
#            Summary          #
###############################
SummaryTableN <- rbind(cbind(FailTime10Output.NoCensor[1],FailTime10Output.NoCensor[2],FailTime10Output.NoCensor[3],FailTime10Output.NoCensor[4],FailTime10Output.NoCensor[5]),
cbind(FailTime20Output.NoCensor[1],FailTime20Output.NoCensor[2],FailTime20Output.NoCensor[3],FailTime20Output.NoCensor[4],FailTime20Output.NoCensor[5]),
cbind(FailTime30Output.NoCensor[1],FailTime30Output.NoCensor[2],FailTime30Output.NoCensor[3],FailTime30Output.NoCensor[4],FailTime30Output.NoCensor[5]),
cbind(FailTime40Output.NoCensor[1],FailTime40Output.NoCensor[2],FailTime40Output.NoCensor[3],FailTime40Output.NoCensor[4],FailTime40Output.NoCensor[5]),
cbind(MedFailTimeOutput.NoCensor[1],MedFailTimeOutput.NoCensor[2],MedFailTimeOutput.NoCensor[3],MedFailTimeOutput.NoCensor[4],MedFailTimeOutput.NoCensor[5]),
cbind("","","","",""),
#cbind(FailTime10Output.Exact[1],FailTime10Output.Exact[2],FailTime10Output.Exact[3],FailTime10Output.Exact[4],FailTime10Output.Exact[5]),
#cbind(FailTime20Output.Exact[1],FailTime20Output.Exact[2],FailTime20Output.Exact[3],FailTime20Output.Exact[4],FailTime20Output.Exact[5]),
#cbind(FailTime30Output.Exact[1],FailTime30Output.Exact[2],FailTime30Output.Exact[3],FailTime30Output.Exact[4],FailTime30Output.Exact[5]),
#cbind(FailTime40Output.Exact[1],FailTime40Output.Exact[2],FailTime40Output.Exact[3],FailTime40Output.Exact[4],FailTime40Output.Exact[5]),
#cbind(MedFailTimeOutput.Exact[1],MedFailTimeOutput.Exact[2],MedFailTimeOutput.Exact[3],MedFailTimeOutput.Exact[4],MedFailTimeOutput.Exact[5]),
#cbind("","","","",""),
#cbind(FailTime10Output.Left[1],FailTime10Output.Left[2],FailTime10Output.Left[3],FailTime10Output.Left[4],FailTime10Output.Left[5]),
#cbind(FailTime20Output.Left[1],FailTime20Output.Left[2],FailTime20Output.Left[3],FailTime20Output.Left[4],FailTime20Output.Left[5]),
#cbind(FailTime30Output.Left[1],FailTime30Output.Left[2],FailTime30Output.Left[3],FailTime30Output.Left[4],FailTime30Output.Left[5]),
#cbind(FailTime40Output.Left[1],FailTime40Output.Left[2],FailTime40Output.Left[3],FailTime40Output.Left[4],FailTime40Output.Left[5]),
#cbind(MedFailTimeOutput.Left[1],MedFailTimeOutput.Left[2],MedFailTimeOutput.Left[3],MedFailTimeOutput.Left[4],MedFailTimeOutput.Left[5]),
#cbind("","","","",""),
#cbind(FailTime10Output.Interval[1],FailTime10Output.Interval[2],FailTime10Output.Interval[3],FailTime10Output.Interval[4],FailTime10Output.Interval[5]),
#cbind(FailTime20Output.Interval[1],FailTime20Output.Interval[2],FailTime20Output.Interval[3],FailTime20Output.Interval[4],FailTime20Output.Interval[5]),
#cbind(FailTime30Output.Interval[1],FailTime30Output.Interval[2],FailTime30Output.Interval[3],FailTime30Output.Interval[4],FailTime30Output.Interval[5]),
#cbind(FailTime40Output.Interval[1],FailTime40Output.Interval[2],FailTime40Output.Interval[3],FailTime40Output.Interval[4],FailTime40Output.Interval[5]),
#cbind(MedFailTimeOutput.Interval[1],MedFailTimeOutput.Interval[2],MedFailTimeOutput.Interval[3],MedFailTimeOutput.Interval[4],MedFailTimeOutput.Interval[5]),
#cbind("","","","",""),
#cbind("","","","",""),
cbind(FailTime10bootOutput.NoCensor[1],FailTime10bootOutput.NoCensor[2],FailTime10bootOutput.NoCensor[3],FailTime10bootOutput.NoCensor[4],FailTime10bootOutput.NoCensor[5]),
cbind(FailTime20bootOutput.NoCensor[1],FailTime20bootOutput.NoCensor[2],FailTime20bootOutput.NoCensor[3],FailTime20bootOutput.NoCensor[4],FailTime20bootOutput.NoCensor[5]),
cbind(FailTime30bootOutput.NoCensor[1],FailTime30bootOutput.NoCensor[2],FailTime30bootOutput.NoCensor[3],FailTime30bootOutput.NoCensor[4],FailTime30bootOutput.NoCensor[5]),
cbind(FailTime40bootOutput.NoCensor[1],FailTime40bootOutput.NoCensor[2],FailTime40bootOutput.NoCensor[3],FailTime40bootOutput.NoCensor[4],FailTime40bootOutput.NoCensor[5]),
cbind(MedFailTimebootOutput.NoCensor[1],MedFailTimebootOutput.NoCensor[2],MedFailTimebootOutput.NoCensor[3],MedFailTimebootOutput.NoCensor[4],MedFailTimebootOutput.NoCensor[5])
#cbind("","","","",""),
#cbind(FailTime10bootOutput.Exact[1],FailTime10bootOutput.Exact[2],FailTime10bootOutput.Exact[3],FailTime10bootOutput.Exact[4],FailTime10bootOutput.Exact[5]),
#cbind(FailTime20bootOutput.Exact[1],FailTime20bootOutput.Exact[2],FailTime20bootOutput.Exact[3],FailTime20bootOutput.Exact[4],FailTime20bootOutput.Exact[5]),
#cbind(FailTime30bootOutput.Exact[1],FailTime30bootOutput.Exact[2],FailTime30bootOutput.Exact[3],FailTime30bootOutput.Exact[4],FailTime30bootOutput.Exact[5]),
#cbind(FailTime40bootOutput.Exact[1],FailTime40bootOutput.Exact[2],FailTime40bootOutput.Exact[3],FailTime40bootOutput.Exact[4],FailTime40bootOutput.Exact[5]),
#cbind(MedFailTimebootOutput.Exact[1],MedFailTimebootOutput.Exact[2],MedFailTimebootOutput.Exact[3],MedFailTimebootOutput.Exact[4],MedFailTimebootOutput.Exact[5]),
#cbind("","","","",""),
#cbind(FailTime10bootOutput.Left[1],FailTime10bootOutput.Left[2],FailTime10bootOutput.Left[3],FailTime10bootOutput.Left[4],FailTime10bootOutput.Left[5]),
#cbind(FailTime20bootOutput.Left[1],FailTime20bootOutput.Left[2],FailTime20bootOutput.Left[3],FailTime20bootOutput.Left[4],FailTime20bootOutput.Left[5]),
#cbind(FailTime30bootOutput.Left[1],FailTime30bootOutput.Left[2],FailTime30bootOutput.Left[3],FailTime30bootOutput.Left[4],FailTime30bootOutput.Left[5]),
#cbind(FailTime40bootOutput.Left[1],FailTime40bootOutput.Left[2],FailTime40bootOutput.Left[3],FailTime40bootOutput.Left[4],FailTime40bootOutput.Left[5]),
#cbind(MedFailTimebootOutput.Left[1],MedFailTimebootOutput.Left[2],MedFailTimebootOutput.Left[3],MedFailTimebootOutput.Left[4],MedFailTimebootOutput.Left[5]),
#cbind("","","","",""),
#cbind(FailTime10bootOutput.Interval[1],FailTime10bootOutput.Interval[2],FailTime10bootOutput.Interval[3],FailTime10bootOutput.Interval[4],FailTime10bootOutput.Interval[5]),
#cbind(FailTime20bootOutput.Interval[1],FailTime20bootOutput.Interval[2],FailTime20bootOutput.Interval[3],FailTime20bootOutput.Interval[4],FailTime20bootOutput.Interval[5]),
#cbind(FailTime30bootOutput.Interval[1],FailTime30bootOutput.Interval[2],FailTime30bootOutput.Interval[3],FailTime30bootOutput.Interval[4],FailTime30bootOutput.Interval[5]),
#cbind(FailTime40bootOutput.Interval[1],FailTime40bootOutput.Interval[2],FailTime40bootOutput.Interval[3],FailTime40bootOutput.Interval[4],FailTime40bootOutput.Interval[5]),
#cbind(MedFailTimebootOutput.Interval[1],MedFailTimebootOutput.Interval[2],MedFailTimebootOutput.Interval[3],MedFailTimebootOutput.Interval[4],MedFailTimebootOutput.Interval[5])
)
write.csv(SummaryTableN, file = "C:/Users/HEMI/Documents/GitHub/Programming_Demo/R/SummaryTableN.csv")


}


sum(Weibull.FailTime10SE.NoCensor.Successcount)
sum(Weibull.FailTime10CIRange.NoCensor.Successcount)
sum(Weibull.FailTime20SE.NoCensor.Successcount)
sum(Weibull.FailTime20CIRange.NoCensor.Successcount)
sum(Weibull.FailTime30SE.NoCensor.Successcount)
sum(Weibull.FailTime30CIRange.NoCensor.Successcount)
sum(Weibull.FailTime40SE.NoCensor.Successcount)
sum(Weibull.FailTime40CIRange.NoCensor.Successcount)
sum(Weibull.MedFailTimeSE.NoCensor.Successcount)
sum(Weibull.MedFailTimeCIRange.NoCensor.Successcount)

#sum(Weibull.FailTime10SE.Exact.Successcount)
#sum(Weibull.FailTime10CIRange.Exact.Successcount)
#sum(Weibull.FailTime20SE.Exact.Successcount)
#sum(Weibull.FailTime20CIRange.Exact.Successcount)
#sum(Weibull.FailTime30SE.Exact.Successcount)
#sum(Weibull.FailTime30CIRange.Exact.Successcount)
#sum(Weibull.FailTime40SE.Exact.Successcount)
#sum(Weibull.FailTime40CIRange.Exact.Successcount)
#sum(Weibull.MedFailTimeSE.Exact.Successcount)
#sum(Weibull.MedFailTimeCIRange.Exact.Successcount)

#sum(Weibull.FailTime10SE.Left.Successcount)
#sum(Weibull.FailTime10CIRange.Left.Successcount)
#sum(Weibull.FailTime20SE.Left.Successcount)
#sum(Weibull.FailTime20CIRange.Left.Successcount)
#sum(Weibull.FailTime30SE.Left.Successcount)
#sum(Weibull.FailTime30CIRange.Left.Successcount)
#sum(Weibull.FailTime40SE.Left.Successcount)
#sum(Weibull.FailTime40CIRange.Left.Successcount)
#sum(Weibull.MedFailTimeSE.Left.Successcount)
#sum(Weibull.MedFailTimeCIRange.Left.Successcount)

#sum(Weibull.FailTime10SE.Interval.Successcount)
#sum(Weibull.FailTime10CIRange.Interval.Successcount)
#sum(Weibull.FailTime20SE.Interval.Successcount)
#sum(Weibull.FailTime20CIRange.Interval.Successcount)
#sum(Weibull.FailTime30SE.Interval.Successcount)
#sum(Weibull.FailTime30CIRange.Interval.Successcount)
#sum(Weibull.FailTime40SE.Interval.Successcount)
#sum(Weibull.FailTime40CIRange.Interval.Successcount)
#sum(Weibull.MedFailTimeSE.Interval.Successcount)
#sum(Weibull.MedFailTimeCIRange.Interval.Successcount)