#-----------------------------------------------------------------------------#
#               RAKING WITH IMPUTED DATA                    #
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
library(miceadds)
library(survey)
library(mitools)

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

###--------------------------------raking------------------------------------###

#create a survey design object
des<-svydesign(id=~1, data=imputationList(micecompletedata))

#copy that object into an another to apply weights
small.svy.rake<-des


# assign the new weighted proportions for each levels of the factor
#here we take household density again
hhdenscat<-c(0.721,0.279)


#create a dataframe with the new weighted distribution
hhdenscat.dist <- data.frame(hhdenscat=  c( "<2 pp bedroom",">=2 pp bedroom"),
                             Freq = nrow(des) * c(0.721,0.279))

#loop through the m imputed designs applying the `rake` function to each of them
# with the new weighted distribution
small.svy.rake$designs <- 
  lapply( 
    des$designs , 
    rake ,
    sample.margins = list(~hhdenscat),
    population.margins = list(hhdenscat.dist)
  )

#compute new estimated seroprevalence with svymean() function
MIcombine(with(small.svy.rake,svymean(~sero_result)))
