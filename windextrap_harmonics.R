########################################################################
# windextrap_harmonics.R                              
# Created by M.Alifa  
# 1st version: initial draft for publication, October 20, 2020 
########################################################################
#   This function fits a nonlinear, temporally-varying model for vertical extrapolation of hourly wind speed for the purpose 
#   of wind energy applications, and performs wind speed extrapolation in the presence of the necessary data
#   More details about the model can be found in Crippa et al, 2020 [insert DOI once we have it]
#
# Inputs: 
#
# 1."train" is the dataframe with data for model training, with column names: 
#   "hour" for the hour of day for each observation
#   "v1"  for time series of near ground wind speeds 
#   "h1" for height above ground of each v1 observation
#   "v2"  for time series of wind speeds aloft 
#   "h2" for height above ground of each v2 observation
#
#2. "test" (optional), a dataframe with the same column names as "train" ("v2" is not necessary)
#   if specified, the function will return extrapolated wind speeds and their 95% Confidence interval for the testing set data
#
# Output
# A list with the following variables: 
# 1. "model": class 'nls', harmonics model for wind speed extrapolation
# 2. "res_model": class 'nls', harmonics model for the error variance 
# Additional output when "test" is specified:
# 3. "v2pred": vector of extrapolated wind speeds 
# 4. "v2_lower95": vector of lower 95% CI for the extrapolated wind speed
# 5. "v2_upper95": vector of upper 95% CI for the extrapolated wind speed
#
#
# Function requires the following R libraries: nlstools
#
# Other potential changes: Incorporate flexible K, flexible CI,  
# allow user input of model initial values, allow user input of additional covariate. 
########################################################################

windextrap_harmonics=function(train, test=NULL){
require(nlstools)
  # NLS model starting values 
  #(ballparked by looking at harmonics of a couple of sites, there might be a better way of estimating these?)
  a.0 <- 0.1
  b1.0 <-  -0.01
  b2.0 <- 0.1
  b3.0 <- -0.001
  b4.0 <- -0.01
  b5.0 <-  0.01
  b6.0 <- -0.01
  b7.0 <- -0.01
  b8.0 <- 0.01
  b9.0 <-  -0.01
  b10.0 <- 0.01
#Starting values list:  
start1 <- list(a=a.0, b1=b1.0, b2=b2.0, b3=b3.0, b4=b4.0, b5=b5.0, b6=b6.0, b7=b7.0, b8=b8.0, b9=b9.0, b10=b10.0)

# NLS model for errors -> Inital values taken from sample fit of K=5 error harmonics
int.0 =0
r1.0= -0.001064
r2.0=  1.768218
r3.0=  0.360004
r4.0= -0.372524
r5.0=  -0.142459
r6.0= -0.426044
r7.0=  -0.206310
r8.0=  0.443649
r9.0=  0.175627
r10.0= -0.191680
#Create starting values list for error variance
start2=list(int=int.0, r1=r1.0, r2=r2.0, r3=r3.0, r4=r4.0,r5=r5.0, r6=r6.0,r7=r7.0, r8=r8.0,r9=r9.0, r10=r10.0)

# Catch errors for incorrect input of "train"
datanames=c("v1","v2","h1","h2","hour")
if(class(train)!= "data.frame" | sum(datanames %in% colnames(train))!=5){
  stop("Error: input for 'train' must be type 'data.frame' with column names 'v1','v2','h1','h2','hour'")
  }else{
## Compute weighted least squares harmonics model on the wind speeds ##
hVAR <- tapply(train$v2,train$hour, var) #Hourly variance for use in weighted least squares model
for(h in unique(train$hour)){ # Assign weights to corresponding hours in train dataframe
train$hVAR[(train$hour==h)]=hVAR[which(as.numeric(names(hVAR))==h)]
}
#Catch errors in fitting of the nls function
error1=try(nls(v2 ~ v1  * (h2/h1)^(a)*  (h2/h1)^(b1*sin(2*pi*hour/24) + b2*cos(2*pi*hour/24) + b3*sin(2*pi*hour*2/24) + b4*cos(2*pi*hour*2/24)
                                                 + b5*sin(2*pi*hour*3/24) + b6*cos(2*pi*hour*3/24) + b7*sin(2*pi*hour*4/24) + b8*cos(2*pi*hour*4/24)
                                                 + b9*sin(2*pi*hour*5/24) + b10*cos(2*pi*hour*5/24) )
               , data=train, start = start1, weights=1/train$hVAR
               , control = nls.control( minFactor = 1/2048, maxiter=1000)), silent=TRUE)
if(class(error1)=="try-error"){  #Error message if windspeed model fit has errors
  stop(paste("Error in model fit: ",geterrmessage()))
  }else{ # If nls function works, continue to fit model
  model <- nls(v2 ~ v1  * (h2/h1)^(a)* (h2/h1)^(b1*sin(2*pi*hour/24) + b2*cos(2*pi*hour/24) + b3*sin(2*pi*hour*2/24) + b4*cos(2*pi*hour*2/24)
                                                + b5*sin(2*pi*hour*3/24) + b6*cos(2*pi*hour*3/24) + b7*sin(2*pi*hour*4/24) + b8*cos(2*pi*hour*4/24)
                                                + b9*sin(2*pi*hour*5/24) + b10*cos(2*pi*hour*5/24) )
               , data=train, start = start1,weights=1/train$hVAR
               , control = nls.control( minFactor = 1/2048, maxiter=1000))

  out=list(model=model) # Initialize output list with the wind extrapolation model
  
  ## Fitting model for error variance ##
  train$residuals=residuals(model)
  VAR_RES <- tapply(train$residuals,train$hour, var) #Calculate residual variances for each hour
  hour=as.numeric(names(VAR_RES)) # Get corresponding hours for residual model input
  
  # Catch errors in fitting of the nls function on residuals
  error2=try(nls(log(VAR_RES) ~ int+ r1*sin(2*pi*hour/24) + r2*cos(2*pi*hour/24) + r3*sin(2*pi*hour*2/24) + r4*cos(2*pi*hour*2/24)
                 + r5*sin(2*pi*hour*3/24) + r6*cos(2*pi*hour*3/24) + r7*sin(2*pi*hour*4/24) + r8*cos(2*pi*hour*4/24)
                 + r9*sin(2*pi*hour*5/24) + r10*cos(2*pi*hour*5/24)
                 , start = start2,
                 control = nls.control( minFactor = 1/2048, maxiter=1000)), silent=TRUE)
  
  if(class(error2)=="try-error"){ # Error message for problems in residual model fit
    print(paste("Error in residual model fit: ",geterrmessage()), "function will not output confidence intervals either")
    }else{ # If nls function works, continue to fit model
    resvar_model <- nls(log(VAR_RES) ~ int+ r1*sin(2*pi*hour/24) + r2*cos(2*pi*hour/24) + r3*sin(2*pi*hour*2/24) + r4*cos(2*pi*hour*2/24)
                        + r5*sin(2*pi*hour*3/24) + r6*cos(2*pi*hour*3/24) + r7*sin(2*pi*hour*4/24) + r8*cos(2*pi*hour*4/24)
                        + r9*sin(2*pi*hour*5/24) + r10*cos(2*pi*hour*5/24)
                        , start = start2,
                        control = nls.control( minFactor = 1/2048, maxiter=1000))
    
    out$res_model=resvar_model # Add residual model to output list
    }# end fit of residual model
  
  ## Extrapolating wind speeds ##
    if(!is.null(test)){  # "test" has to be specified for extrapolation routine to happen
      # Catch errors for incorrect input of "test"
      if(class(test)!= "data.frame" | sum((datanames[-2]) %in% colnames(test))!=4){
        print("Error: input for 'test' must be type 'data.frame' with column names 'v1','h1','h2','hour', function will not output wind speed predictions")
        }else{
        out$v2pred=predict(model, newdata=test) # Add predicted wind speeds to the output list
        
        if(class(error2)!="try-error"){ # If residual variance model was fit, calculate CI
        # Calculating confidence interval #
        VAR_RES_MOD=exp(predict(resvar_model)) # Modeled residual variance for the available hours
        for(h in unique(test$hour)){
          #Compute standard deviation and assign to corresponding hours in "test" dataframe
          test$res_sd[(test$hour==h)]=sqrt(VAR_RES_MOD[which(as.numeric(names(VAR_RES))==h)])
        }
        # Use hourly standard deviation to calculate confidence interval
        test$v2pred=predict(model, newdata=test) # Add predicted wind speeds to test dataframe
        lwr=test$v2pred-1.96*test$res_sd #Lower interval
        uppr=test$v2pred+1.96*test$res_sd #Upper interval
        out$v2_lower95=lwr ;  out$v2_upper95=uppr   # Add CI to output list   
        } # end of CI calculation 
      }# end of predictions using "test"
    } #!is.null (test)
  
  return(out) #Return output list
  
  }# end of wind speed model fit
 } # all data in "train" is present
}# end of function
