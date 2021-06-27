library(dplyr)
library(lubridate)
library(readr)
library(deSolve)
library(rstudioapi)
library(zoo)
library(tidyr)

directory <- dirname(rstudioapi::getSourceEditorContext()$path) 



latent_period <- 5
interval=30
infectious_period <- 9

#data conversion
cases=read.csv(file.path(directory, "data/covid19-download.csv"))
#cases <- read.csv("covid19-download.csv",header=T)
mydata = cases %>% filter(pruid==35) %>%
        select(prname,date,numtoday,numdeathstoday,numrecoveredtoday,numactive,numtestedtoday,numteststoday,numtests,numrecover) %>%
        rowwise() %>% mutate(numtestsmax=max(numtestedtoday,numteststoday,na.rm=T))%>%
        select(-numtestedtoday,-numteststoday) %>%
        drop_na() %>% mutate(date=as_date(date)) %>%
        mutate(conftotal=cumsum(numtoday)) %>%
        mutate(S=max(numtests)) %>%
        mutate(E=rollapplyr(numtoday,width=latent_period,FUN=sum,partial=T)) %>%
        mutate(I=numactive-E ,R=numrecover) %>% drop_na(conftotal)%>%
        complete(date = full_seq(date, period = 1))
# mydata=cases[cases$pruid==35,]
# mydata=mydata[,which(names(mydata) %in% c("prname","date","numtoday","numdeathstoday","numrecoveredtoday","numactive","numtestedtoday","numteststoday"))]
# mydata2=mydata[,which(names(mydata) %in% c("numtestedtoday","numteststoday"))]
# mydata$numtests=apply(mydata2,1,max,na.rm=T)
# mydata[mydata$numtests==-Inf,]$numtests=NA
# mydata=mydata[,-which(names(mydata) %in% c("numtestedtoday","numteststoday"))]
# mydata=na.omit(mydata)
# mydata$date=as_date(mydata$date)
start=mydata$date[1]
mydata$time=mydata$date-start+1
#n=2.3 * 1000000/1000
mydata$S=mydata$S-mydata$E-mydata$I-mydata$numrecover
# mydata$E=mydata$numtoday+within 5 days
# mydata$I=mydata$numactive-mydata$numtoday-within 5 days
# mydata$R=mydata$numrecoveredtoday
mydata$N=mydata$S+mydata$E+mydata$I+mydata$R


#function to calculate beta

sir.calc <- function(date, N=14570000,t=1){
  Inew <- mydata[mydata$date == date, c("numtoday")]  # number of new infections since last sample
  I <- mydata[mydata$date == date, c("numactive")] # number of infected individuals
  R<- mydata[mydata$date == date, c("R")] # number of recovered 
  D<- mydata[mydata$date == date, c("numdeathstoday")] # number of deaths
  S <- mydata[mydata$date == date, c("S")]
  beta <- (-1/t) * log(1-Inew*((1/I)+(1/S)))
}

mydata$beta=mapply(sir.calc,mydata$date,mydata$N)
mydata$beta


##seir function
seir <- function(beta, gamma, delta, S, E, I, R,i) {
  parameter_list <- c(beta, gamma, delta)
  W <- S
  X <- E
  Y <- I
  Z <- R
  N<- S+E+I+R
  
  seir_model = function (current_timepoint, state_values, parameters)
  {
    # create state variables (local variables)
    S = state_values [1]        # susceptibles
    E = state_values [2]        # exposed
    I = state_values [3]        # infectious
    R = state_values [4]        # recovered
    
    with ( 
      as.list (parameters),     # variable names within parameters can be used 
      {
        # compute derivatives
        dS = (-beta * S * I) 
        dE = (beta * S * I) - (delta * E) 
        dI = (delta * E) - (gamma * I) 
        dR = (gamma * I) 
        
        
        # combine results
        results = c (dS, dE, dI, dR)
        list (results)
      }
    )
  }
  
  initial_values = c (S = W/N, E = X/N, I = Y/N, R = Z/N)
  timepoints = seq (i, i+interval-1, by=1)
  output = lsoda(initial_values, timepoints, seir_model, parameter_list)
  return(output)
}

seir2 <- function(beta, gamma, delta, S, E, I, R,D,i) {
  parameter_list <- c(beta, gamma, delta)
  W <- S
  X <- E
  Y <- I
  Z <- R
  N<- S+E+I+R
  D= D/N
  
  seir_model2 = function (current_timepoint, state_values, parameters)
  {
    # create state variables (local variables)
    S = state_values [1]        # susceptibles
    E = state_values [2]        # exposed
    I = state_values [3]        # infectious
    R = state_values [4]        # recovered
    
    with ( 
      as.list (parameters),     # variable names within parameters can be used 
      {
        # compute derivatives
        # dS = (-beta * S * I) 
        # dE = (beta * S * I) - (delta * E) 
        # dI = (delta * E) - (gamma * I) 
        # dR = (gamma * I) 
        dS = (-beta * S * I) - (D * S)
        dE = (beta * S * I) - (delta * E) - (D * I)
        dI = (delta * E) - (gamma * I) - (D * E)
        dR = (gamma * I) - (D * R)
        
        # combine results
        results = c (dS, dE, dI, dR)
        list (results)
      }
    )
  }
  
  initial_values = c (S = W/N, E = X/N, I = Y/N, R = Z/N)
  timepoints = seq (i, i+interval-1, by=1)
  output = lsoda(initial_values, timepoints, seir_model2, parameter_list)
  return(output)
}


#parameter setting
#N <- 2.3 * 1000000

t=7

gamma <- 1/infectious_period
delta <- 1/latent_period
#beta <- (-1/t) * log(1-Inew*((1/I)+(1/S)))
#beta=0.111

#output_test <- seir(mydata$beta[1], gamma, delta, mydata$S[1], mydata$E[1], mydata$I[1], mydata$R[1],1)
i=seq(1,nrow(mydata),interval)
df=mydata[seq(1,nrow(mydata),interval),]
output=mapply(seir,df$beta,gamma,delta,df$S,df$E,df$I,df$R,i,SIMPLIFY=FALSE)
output=Reduce( rbind.data.frame,output)

output2=mapply(seir2,df$beta,gamma,delta,df$S,df$E,df$I,df$R,df$numdeathstoday,i,SIMPLIFY=FALSE)
output2=Reduce( rbind.data.frame,output2)



#observed plot vs estimated rate
par(mfrow=c(2, 2))
ymax=max(mydata$S/mydata$N,output$S,na.rm=T)
ymin=min(mydata$S/mydata$N,output$S,na.rm=T)
matplot(mydata$time,mydata$S/mydata$N,type = 'l',col='blue',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="rate",main="Susceptible")
par(new = TRUE)
matplot(output$time,output$S,type = 'l',col='red',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="rate")

ymax=max(mydata$E/mydata$N,output$E,na.rm=T)
ymin=min(mydata$E/mydata$N,output$E,na.rm=T)
matplot(mydata$time,mydata$E/mydata$N,type = 'l',col='blue',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="rate",main="Exposed")
par(new = TRUE)
matplot(output$time,output$E,type = 'l',col='red',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="rate")

ymax=max(mydata$I/mydata$N,output$I,na.rm=T)
ymin=min(mydata$I/mydata$N,output$I,na.rm=T)
matplot(mydata$time,mydata$I/mydata$N,type = 'l',col='blue',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="rate",main="Infectious")
par(new = TRUE)
matplot(output$time,output$I,type = 'l',col='red',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="rate")

ymax=max(mydata$R/mydata$N,output$R,na.rm=T)
ymin=min(mydata$R/mydata$N,output$R,na.rm=T)
matplot(mydata$time,mydata$R/mydata$N,type = 'l',col='blue',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="rate",main="Recovered")
par(new = TRUE)
matplot(output$time,output$R,type = 'l',col='red',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="rate")



#observed plot vs estimated num
output=output[1:nrow(mydata),]

par(mfrow=c(2, 2))
ymax=max(mydata$S,output$S*mydata$N,na.rm=T)
ymin=min(mydata$S,output$S*mydata$N,na.rm=T)
#xmax=mydata$date[1]
#xmin=mydata$date[length(mydata$date)]
matplot(mydata$time,mydata$S,type = 'l',col='blue',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="count",main="Susceptible")
par(new = TRUE)
matplot(output$time,output$S*mydata$N,type = 'l',col='red',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="count")

ymax=max(mydata$E,output$E*mydata$N,na.rm=T)
ymin=min(mydata$E,output$E*mydata$N,na.rm=T)
matplot(mydata$time,mydata$E,type = 'l',col='blue',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="count",main="Exposed")
par(new = TRUE)
matplot(output$time,output$E*mydata$N,type = 'l',col='red',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="count")

#interventions=as.Date(c('2020-09-08','2020-09-17','2020-09-25','2020-09-28','2020-10-02','2020-10-10','2020-10-16','2020-11-06','2020-11-07','2020-11-14',
#                       '2020-11-16','2020-11-23','2020-12-14','2020-12-15','2020-12-21','2020-12-26','2021-01-14','2021-01-29','2021-02-10','2021-03-03','2021-03-05','2021-04-08'))
#lockdowns=as.Date(c('2020-12-26','2021-04-08'))
#abline(v=interventions-start+1) 
#abline(v=lockdowns-start+1,col="orange")


ymax=max(mydata$I,output$I*mydata$N,na.rm=T)
ymin=min(mydata$I,output$I*mydata$N,na.rm=T)
matplot(mydata$time,mydata$I,type = 'l',col='blue',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="count",main="Infectious")
par(new = TRUE)
matplot(output$time,output$I*mydata$N,type = 'l',col='red',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="count")

ymax=max(mydata$R,output$R*mydata$N,na.rm=T)
ymin=min(mydata$R,output$R*mydata$N,na.rm=T)
matplot(mydata$time,mydata$R,type = 'l',col='blue',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="count",main="Recovered")
par(new = TRUE)
matplot(output$time,output$R*mydata$N,type = 'l',col='red',ylim=c(ymin,ymax),xlim=c(0,410),xlab="time",ylab="count")


interventions=as.Date(c('2020-09-08','2020-09-17','2020-09-25','2020-09-28','2020-10-02','2020-10-10','2020-10-16','2020-11-06','2020-11-07','2020-11-14',
                        '2020-11-16','2020-11-23','2020-12-14','2020-12-15','2020-12-21','2020-12-26','2021-01-14','2021-01-29','2021-02-10','2021-03-03','2021-03-05','2021-04-08'))
lockdowns=as.Date(c('2020-12-26','2021-04-08'))
abline(v=interventions-start+1) 
abline(v=lockdowns-start+1,col="orange")

#rate ggplot
names(output)=c("time","Shat","Ehat","Ihat","Rhat")
mydata$time=as.numeric(mydata$time)
df=inner_join(mydata,output,by="time")
end=mydata$date[length(mydata$date)]
p1=ggplot(df,aes(x=date))+geom_line(aes(y=S/N,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Susceptible (rate)", labels = scales::comma)+
  geom_line(mapping=aes(y=Shat,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank())
p2=ggplot(df,aes(x=date))+geom_line(aes(y=E/N,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Exposed (rate)", labels = scales::comma)+
  geom_line(mapping=aes(y=Ehat,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank())
p3=ggplot(df,aes(x=date))+geom_line(aes(y=I/N,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Infectious (rate)", labels = scales::comma)+
  geom_line(mapping=aes(y=Ihat,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank())
p4=ggplot(df,aes(x=date))+geom_line(aes(y=R/N,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Recovered (rate)", labels = scales::comma)+
  geom_line(mapping=aes(y=Rhat,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank())

library(ggpubr)
ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

  
#count ggplot
p1=ggplot(df,aes(x=date))+geom_line(aes(y=S,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Susceptible (count)", labels = scales::comma)+
  geom_line(mapping=aes(y=Shat*N,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank()) 
p2=ggplot(df,aes(x=date))+geom_line(aes(y=E,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Exposed (count)", labels = scales::comma)+
  geom_line(mapping=aes(y=Ehat*N,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank()) 
p3=ggplot(df,aes(x=date))+geom_line(aes(y=I,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Infectious (count)", labels = scales::comma)+
  geom_line(mapping=aes(y=Ihat*N,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank()) 
p4=ggplot(df,aes(x=date))+geom_line(aes(y=R,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Recovered (count)", labels = scales::comma)+
  geom_line(mapping=aes(y=Rhat*N,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank()) 

ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

#function to calculated different measures of error
rmse <- function(pred, obs){
  mse <- mean((pred - obs)^2, na.rm = TRUE)
  rmse <- sqrt(mse)
  rrmse <- rmse/sd(obs, na.rm = TRUE)
  print(c(mse, rmse, rrmse))
}

# calculate values for S, E, I, R 
rmse.s <- rmse(predicted$S, observed$S) #1.498821e+08 1.224263e+04 8.365103e-02
rmse.e <- rmse(predicted$E, observed$E) #1.378359e+07 3.712626e+03 6.375511e-01
rmse.i <- rmse(predicted$I, observed$I) #6.751442e+06 2.598354e+03 4.812551e-01
rmse.r <- rmse(predicted$R, observed$R) #1.014787e+08 1.007366e+04 7.354005e-02



#seir2
# par(mfrow=c(2, 2))
# matplot(mydata$time,mydata$S/mydata$N,type = 'l',col='blue',ylim=c(0,1),xlim=c(0,410),xlab="time",ylab="rate",main="Susceptible")
# par(new = TRUE)
# matplot(output2$time,output2$S,type = 'l',col='red',ylim=c(0,1),xlim=c(0,410),xlab="time",ylab="rate")
# 
# matplot(mydata$time,mydata$E/mydata$N,type = 'l',col='blue',ylim=c(0,1),xlim=c(0,410),xlab="time",ylab="rate",main="Exposed")
# par(new = TRUE)
# matplot(output2$time,output2$E,type = 'l',col='red',ylim=c(0,1),xlim=c(0,410),xlab="time",ylab="rate")
# 
# matplot(mydata$time,mydata$I/mydata$N,type = 'l',col='blue',ylim=c(0,1),xlim=c(0,410),xlab="time",ylab="rate",main="Infectious")
# par(new = TRUE)
# matplot(output2$time,output2$I,type = 'l',col='red',ylim=c(0,1),xlim=c(0,410),xlab="time",ylab="rate")
# 
# matplot(mydata$time,mydata$R/mydata$N,type = 'l',col='blue',ylim=c(0,1),xlim=c(0,410),xlab="time",ylab="rate",main="Recovered")
# par(new = TRUE)
# matplot(output2$time,output2$R,type = 'l',col='red',ylim=c(0,1),xlim=c(0,410),xlab="time",ylab="rate")
# abline(v=as.Date('2020-09-08')-start+1) 
