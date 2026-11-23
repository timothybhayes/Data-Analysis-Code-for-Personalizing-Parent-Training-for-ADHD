##################################################################
# Results Formatting Functions and Result Table Creation 3-12-26 #
##################################################################
rm(list = ls())


dat <- read.csv('moduleDataCorrected.csv')
library("MRTAnalysis")


ciFormat <- function(lower, upper) paste0("[", lower, ", ", upper, "]")
mrtRes <- function(mrtModel, treatName = "Treatment", modNames = NULL){
	res <- summary(mrtModel)$causal_excursion_effect	
	res <- as.data.frame(res)
	rownames(res)[1] <- treatName
	if(!is.null(modNames)){
		rownames(res)[-1] <- paste0(treatName, " x ", modNames)
	}
	res <- res[,c("Estimate", "StdErr", "95% LCL", "95% UCL")]
	res$RR <- exp(res$Estimate)
	res$expLCL <- exp(res[,"95% LCL"])
	res$expUCL <- exp(res[,"95% UCL"])
	
	res[,c("Estimate", "RR")] <- round(res[,c("Estimate", "RR")], 3)
	res[,c("StdErr", "95% LCL", "95% UCL", "expLCL", "expUCL")] <- round(res[,c("StdErr", "95% LCL", "95% UCL")], 2)
	
	res$CI <- res$expCI <-  NA
	for(i in 1:nrow(res)){
		res$CI[i] <- ciFormat(res[i,"95% LCL"], res[i,"95% UCL"])
		res$expCI[i] <- ciFormat(res[i,"expLCL"], res[i,"expUCL"])
	}
	
	res <- res[,c("Estimate", "StdErr", "CI", "RR", "expCI")]
	colnames(res) <- c("b", "SE", "95% CI", "RR", "RR 95% CI")
	return(res)
}


#***********************************#
# Naming Conventions Moving Forward #
#***********************************#

#ct_ = control vs treatment (prefix)

#_p or _n = positive v other or negative v other respectively

######################################################
# 1. Full Sample (Not-As-Exploratory of =p) Analyses #
######################################################

#=============================#
# A. Causal Excursion Effects #
#=============================#

# Note to self: center various moderators so that 0 is meaningful

#Create probability column with .5 for intervention vs. no intervention
dat$prob <- .5


ct_excursion_p <- emee(
    data = na.omit(dat[,c("userid", "positive_out", "Intervention", "prob")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~1,
    control_formula = ~1
    )

summary(ct_excursion_p)

ct_excursion_n <- emee(
    data = na.omit(dat[,c("userid", "negative_out", "Intervention", "prob")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~1,
    control_formula = ~1
    )

summary(ct_excursion_n)

#========================#
# B. Covariate By Itself #
#========================#

#centering age for these analyses ... do not save this data file to overwrite any existing one.
dat$age <- scale(dat$age, scale = FALSE)

ct_covariate_p <- emee(
    data = na.omit(dat[,c("userid", "positive_out", "Intervention", "age", "prob")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~1,
    control_formula = ~ age
    )

summary(ct_covariate_p)

ct_covariate_n <- emee(
    data = na.omit(dat[,c("userid", "negative_out", "Intervention", "age", "prob")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~1,
    control_formula = ~ age
    )

summary(ct_covariate_n)

#=======================#
# C. Moderator Analyses #
#=======================#


#***************#
# Module Number #
#***************#

dat$modulenumc <- scale(dat$modulenum, scale = FALSE) 

ct_modulenum_p <- emee(
    data = na.omit(dat[,c("userid", "positive_out", "Intervention", "age", "prob", "modulenumc")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~modulenumc,
    control_formula = ~modulenumc + age
    )

summary(ct_modulenum_p)

ct_modulenum_n <- emee(
    data = na.omit(dat[,c("userid", "negative_out", "Intervention", "age", "prob", "modulenumc")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~modulenumc,
    control_formula = ~modulenumc + age
    )

summary(ct_modulenum_n)

#***************#
# Lagged Affect #
#***************#

ct_laggedaffect_p <- emee(
    data = na.omit(dat[,c("userid", "positive_out", "Intervention", "age", "prob", "positive", "neutral")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~positive + neutral,
    control_formula = ~positive + neutral + age
    )

summary(ct_laggedaffect_p)

ct_laggedaffect_n <- emee(
    data = na.omit(dat[,c("userid", "negative_out", "Intervention", "age", "prob", "positive", "neutral")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~positive + neutral,
    control_formula = ~positive + neutral + age
    )

summary(ct_laggedaffect_n)

#*****#
# Age #
#*****#

ct_agemod_p <- emee(
    data = na.omit(dat[,c("userid", "positive_out", "Intervention", "age", "prob")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~age,
    control_formula = ~ age
    )

summary(ct_agemod_p)

ct_agemod_n <- emee(
    data = na.omit(dat[,c("userid", "negative_out", "Intervention", "age", "prob")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~age,
    control_formula = ~ age
    )

summary(ct_agemod_n)

#*****#
# ODD #
#*****#

#Oh, Dee Dee ...

dat$ODDsevc <- scale(dat$ODDsev, scale = FALSE)


ct_odd_p <- emee(
    data = na.omit(dat[,c("userid", "positive_out", "Intervention", "age", "prob", "ODDsevc")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~ODDsevc,
    control_formula = ~ODDsevc + age
    )

summary(ct_odd_p)

ct_odd_n <- emee(
    data = na.omit(dat[,c("userid", "negative_out", "Intervention", "age", "prob", "ODDsevc")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~ODDsevc,
    control_formula = ~ ODDsevc + age
    )

summary(ct_odd_n)

#******#
# phq9 #
#******#

dat$phq9_scorec <- scale(dat$phq9_score, scale = FALSE)



ct_phq_p <- emee(
    data = na.omit(dat[,c("userid", "positive_out", "Intervention", "age", "prob", "phq9_scorec")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~phq9_scorec,
    control_formula = ~phq9_scorec + age
    )

summary(ct_phq_p)

ct_phq_n <- emee(
    data = na.omit(dat[,c("userid", "negative_out", "Intervention", "age", "prob", "phq9_scorec")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "Intervention",
    rand_prob = "prob",
    moderator_formula = ~phq9_scorec,
    control_formula = ~ phq9_scorec + age
    )

summary(ct_phq_n)

blankRow <- rep(NA, 10)

ct_resTable <- rbind(model1 = blankRow, cbind(mrtRes(ct_excursion_p), mrtRes(ct_excursion_n)), model2 = blankRow, cbind(mrtRes(ct_covariate_p), mrtRes(ct_covariate_n)), model3 = blankRow, cbind(mrtRes(ct_modulenum_p, modNames = "Modules"), mrtRes(ct_modulenum_n, modNames = "Modules")), model4 = blankRow, cbind(mrtRes(ct_laggedaffect_p, modNames = c("Positive", "Neutral")), mrtRes(ct_laggedaffect_n, modNames = c("Positive", "Neutral"))), model5 = blankRow, cbind(mrtRes(ct_agemod_p, modNames = "age"), mrtRes(ct_agemod_n, modNames = "age")), model6 = blankRow, cbind(mrtRes(ct_odd_p, modNames = "ODD"), mrtRes(ct_odd_n, modNames = "ODD")), model7 = blankRow, cbind(mrtRes(ct_phq_p, modNames = "Depression"), mrtRes(ct_phq_n, modNames = "Depression")))

###########################
# 2. Exploratory Analyses #
###########################

dat$StrategyDummy <- ifelse(dat$Intervention.Type == "Parenting Strategy", 1, 0)

datSub <- dat[dat$Intervention.Type != "No Intervention" & dat$Intervention.Type != "",]

#=============================#
# A. Causal Excursion Effects #
#=============================#

# Note to self: center various moderators so that 0 is meaningful


tt_excursion_p <- emee(
    data = na.omit(datSub[,c("userid", "positive_out", "StrategyDummy", "prob")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~1,
    control_formula = ~1
    )

summary(tt_excursion_p)

tt_excursion_n <- emee(
    data = na.omit(datSub[,c("userid", "negative_out", "StrategyDummy", "prob")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~1,
    control_formula = ~1
    )

summary(tt_excursion_n)

#========================#
# B. Covariate By Itself #
#========================#

#centering age for these analyses ... do not save this data file to overwrite any existing one.
dat$age <- scale(dat$age, scale = FALSE)

tt_covariate_p <- emee(
    data = na.omit(datSub[,c("userid", "positive_out", "StrategyDummy", "age", "prob")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~1,
    control_formula = ~ age
    )

summary(tt_covariate_p)

tt_covariate_n <- emee(
    data = na.omit(datSub[,c("userid", "negative_out", "StrategyDummy", "age", "prob")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~1,
    control_formula = ~ age
    )

summary(tt_covariate_n)

#=======================#
# C. Moderator Analyses #
#=======================#


#***************#
# Module Number #
#***************#



tt_modulenum_p <- emee(
    data = na.omit(datSub[,c("userid", "positive_out", "StrategyDummy", "age", "prob", "modulenumc")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~modulenumc,
    control_formula = ~modulenumc + age
    )

summary(tt_modulenum_p)

tt_modulenum_n <- emee(
    data = na.omit(datSub[,c("userid", "negative_out", "StrategyDummy", "age", "prob", "modulenumc")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~modulenumc,
    control_formula = ~modulenumc + age
    )

summary(tt_modulenum_n)

#***************#
# Lagged Affect #
#***************#

tt_laggedaffett_p <- emee(
    data = na.omit(datSub[,c("userid", "positive_out", "StrategyDummy", "age", "prob", "positive", "neutral")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~positive + neutral,
    control_formula = ~positive + neutral + age
    )

summary(tt_laggedaffett_p)

tt_laggedaffett_n <- emee(
    data = na.omit(datSub[,c("userid", "negative_out", "StrategyDummy", "age", "prob", "positive", "neutral")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~positive + neutral,
    control_formula = ~positive + neutral + age
    )

summary(tt_laggedaffett_n)

#*****#
# Age #
#*****#

tt_agemod_p <- emee(
    data = na.omit(datSub[,c("userid", "positive_out", "StrategyDummy", "age", "prob")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~age,
    control_formula = ~ age
    )

summary(tt_agemod_p)

tt_agemod_n <- emee(
    data = na.omit(datSub[,c("userid", "negative_out", "StrategyDummy", "age", "prob")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~age,
    control_formula = ~ age
    )

summary(tt_agemod_n)

#*****#
# ODD #
#*****#

tt_odd_p <- emee(
    data = na.omit(datSub[,c("userid", "positive_out", "StrategyDummy", "age", "prob", "ODDsevc")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~ODDsevc,
    control_formula = ~ODDsevc + age
    )

summary(tt_odd_p)

tt_odd_n <- emee(
    data = na.omit(datSub[,c("userid", "negative_out", "StrategyDummy", "age", "prob", "ODDsevc")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~ODDsevc,
    control_formula = ~ ODDsevc + age
    )

summary(tt_odd_n)

#******#
# phq9 #
#******#

tt_phq_p <- emee(
    data = na.omit(datSub[,c("userid", "positive_out", "StrategyDummy", "age", "prob", "phq9_scorec")]),
    id = "userid",
    outcome = "positive_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~phq9_scorec,
    control_formula = ~phq9_scorec + age
    )

summary(tt_phq_p)

tt_phq_n <- emee(
    data = na.omit(datSub[,c("userid", "negative_out", "StrategyDummy", "age", "prob", "phq9_scorec")]),
    id = "userid",
    outcome = "negative_out",
    treatment = "StrategyDummy",
    rand_prob = "prob",
    moderator_formula = ~phq9_scorec,
    control_formula = ~ phq9_scorec + age
    )

summary(tt_phq_n)

#blankRow <- rep(NA, 10)

tt_resTable <- rbind(model1 = blankRow, cbind(mrtRes(tt_excursion_p), mrtRes(tt_excursion_n)), model2 = blankRow, cbind(mrtRes(tt_covariate_p), mrtRes(tt_covariate_n)), model3 = blankRow, cbind(mrtRes(tt_modulenum_p, modNames = "Modules"), mrtRes(tt_modulenum_n, modNames = "Modules")), model4 = blankRow, cbind(mrtRes(tt_laggedaffett_p, modNames = c("Positive", "Neutral")), mrtRes(tt_laggedaffett_n, modNames = c("Positive", "Neutral"))), model5 = blankRow, cbind(mrtRes(tt_agemod_p, modNames = "age"), mrtRes(tt_agemod_n, modNames = "age")), model6 = blankRow, cbind(mrtRes(tt_odd_p, modNames = "ODD"), mrtRes(tt_odd_n, modNames = "ODD")), model7 = blankRow, cbind(mrtRes(tt_phq_p, modNames = "Depression"), mrtRes(tt_phq_n, modNames = "Depression")))

################
# Output Files #
################

write.csv(ct_resTable, file = 'ctResTable.csv', na = "")

write.csv(tt_resTable, file = 'ttResTable.csv', na = "")