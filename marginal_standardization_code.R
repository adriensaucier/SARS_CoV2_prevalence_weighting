#-----------------------------------------------------------------------------#
#               MARGINAL STANDARDIZATION WITH IMPUTED DATA                    #
#-----------------------------------------------------------------------------#

# Adrien Saucier

# This code is part of the study:

#Generalizability of anti-SARS-CoV-2 seroprevalence estimates to the Montréal
#pediatric population: a comparison between two weighting methods

#by Adrien Saucier, Bouchra Nasri, Britt McKinnon, Mabel Carabali,
# Katia Charland, Laura Pierce & Kate Zinszer

###-----------------------------libraries-----------------------------------###

library(tidyverse)
library(mice)
library(marginaleffects)

###-------------------------------data--------------------------------------###

#run chained equation for imputation
options(max.print = 99999)
micenewdata<-mice(encore_data,maxit=0) #dry run

#predictors matrix#
#specify a predictors matrix
pred<-micenewdata$predictorMatrix
pred[,c("hhid")]<- -2 #set household Id as clustering factor in imputation

#rerun chained equation with new prediction matrix
micenewdata<-mice(encore_data,m=20,pred=pred,seed=31)

#extract the m datasets in long and list forms
micelongdata <- complete(micenewdata, "long")
micelongdata$.imp<-factor(micelongdata$.imp)

###----------------------marginal standardization----------------------------###

# In order to weight for categorical variable as a proportion in R with the 
# marginaeffects package we need to set the variable as a numeric (0,1)
# and feed afterward the desired proportion to the predictions function

# modify the variable to be 0s and 1s in the long format
micelongdata<-micelongdata%>%
  mutate(hhdensnum=case_when(
    hhdenscat==">=2 pp bedroom"~1,
    hhdenscat%in%c("<2 pp bedroom")~0
  ))

#return into list format
micecompletedata<-as.list(split(micelongdata,as.factor(micelongdata$.imp)))

#create a function to build m regressions, seting the weighted proportion 
#in grid_type argument is set at "counterfactual" to manipulate the hhdensnum
# variable value. All other variables are let at the observed values in the 
# dataset. As explained by Muller and al. in their article "Estimating predicted 
#probabilities from logistic regression : different methods correspond to 
#different target population" (DOI:10.1093/ije/dyu029):
#"to apply [marignal standardization] in practice after performing a logistic
#regression, the exposure E is set to the (possibly counterfactual) level e for
#everyone in the dataset, and the logistic regression coefficients are used to
#calculate predicted probabilities for everyone at their observed confounder 
#pattern and newly assigned exposure value".


fit_reg <- function(dat) {
  mod <- glm(sero_result~ethnic_minority+as.numeric(hhdensnum)+income,data=dat,family = binomial(link="logit"))
  out <- predictions(mod,datagrid(hhdensnum=0.279,grid_type = "counterfactual"),by=c("hhdensnum"),type="response")
  return(out)
}

#apply the function to the list format data
model_response <- lapply(micecompletedata, fit_reg)

#print the summary, giving the predicted probabilities when household
# density of more than 2 persons per bedroom represent 27.9% of the population
summary(pool(model_response),conf.int=T)

#the same thing can be done with a categorical variable for j levels 
#(more than 2 levels) by creating j-1 auxiliary variables

#create numeric auxiliary variables
micelongdata<-micelongdata%>%
  mutate(bachelors=case_when(
    educat3=="bachelors"~1,
    educat3%in%c("masters or higher","less than bachelors")~0
  ))%>%
  mutate(masters=case_when(
    educat3=="masters or higher"~1,
    educat3%in%c("bachelors","less than bachelors")~0))

#recompose data as a list   
micecompletedata<-as.list(split(micelongdata,as.factor(micelongdata$.imp)))

#create the marginal standardization function
fit_reg <- function(dat) {
  mod <- glm(sero_result~phase3+as.numeric(masters)+as.numeric(bachelors),data=dat,family = binomial(link="logit"))
  out <- predictions(mod,datagrid(bachelors=0.242,masters=0.192,grid_type = "counterfactual"),by=c("bachelors","masters"),type="response")
  return(out)
}

#apply the function to the list data
model_response <- lapply(micecompletedata, fit_reg)

#get the summary
summary(pool(model_response),conf.int=T)
