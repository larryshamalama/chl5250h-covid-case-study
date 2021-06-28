library(dplyr)
library(lubridate)
library(readr)
library(deSolve)
library(rstudioapi)
library(zoo)
library(tidyr)
library(Metrics)

directory <- dirname(rstudioapi::getSourceEditorContext()$path) 



latent_period <- 5
interval=30
infectious_period <- 9

#data conversion
cases=read.csv(file.path(directory, "data/covid19-download.csv"))
cases_v2=read.csv(file.path(directory, "data/covid19-download_v2.csv"))
#cases <- read.csv("covid19-download.csv",header=T)
mydata = cases %>% filter(pruid==35) %>%
        select(prname,date,numtoday,numdeathstoday,numrecoveredtoday,numactive,numtestedtoday,numteststoday,numtests,numrecover) %>%
        rowwise() %>% mutate(numtestsmax=max(numtestedtoday,numteststoday,na.rm=T))%>%
        select(-numtestedtoday,-numteststoday) %>%
        drop_na() %>% mutate(date=as_date(date)) %>%
        mutate(conftotal=cumsum(numtoday)) %>%
        mutate(N=max(numtests)) %>%
        mutate(E=rollapplyr(numtoday,width=latent_period,FUN=sum,partial=T)) %>%
        mutate(I=numactive-E ,R=numrecover) %>% drop_na(conftotal)%>%
        complete(date = full_seq(date, period = 1)) %>% mutate(S=N-E-I-R) 
start=mydata$date[1]
mydata$time=mydata$date-start+1


#function to calculate beta

sir.calc <- function(date, N=14570000,t=1){
  if (which(mydata$date == date)==1){return(NA)}
  else {
  Inew <- mydata[which(mydata$date == date)-1, c("numtoday")]  # number of new infections since last sample
  I <- mydata[mydata$date == date, c("numactive")] # number of infected individuals
  R<- mydata[mydata$date == date, c("R")] # number of recovered 
  D<- mydata[mydata$date == date, c("numdeathstoday")] # number of deaths
  t=as.numeric(mydata[which(mydata$date == date),c("date")]-mydata[which(mydata$date == date)-1,c("date")])
  S <- mydata[mydata$date == date, c("S")]
  beta <- (-1/t) * log(1-Inew*((1/I)+(1/S)))}
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
  scale_y_continuous(name="Susceptible (rate)", labels = scales::comma)+labs(x="Date",size=1.3)+
  geom_line(mapping=aes(y=Shat,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank())
p2=ggplot(df,aes(x=date))+geom_line(aes(y=E/N,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Exposed (rate)", labels = scales::comma)+labs(x="Date",size=1.3)+
  geom_line(mapping=aes(y=Ehat,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank())
p3=ggplot(df,aes(x=date))+geom_line(aes(y=I/N,col="Observed"),size=1.3)+labs(x="Date",size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Infectious (rate)", labels = scales::comma)+labs(x="Date",size=1.3)+
  geom_line(mapping=aes(y=Ihat,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank())
p4=ggplot(df,aes(x=date))+geom_line(aes(y=R/N,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Recovered (rate)", labels = scales::comma)+labs(x="Date",size=1.3)+
  geom_line(mapping=aes(y=Rhat,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank())

library(ggpubr)
ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

  
#count ggplot
p1=ggplot(df,aes(x=date))+geom_line(aes(y=S,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Susceptible (count)", labels = scales::comma)+labs(x="Date",size=1.3)+
  geom_line(mapping=aes(y=Shat*N,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank()) 
p2=ggplot(df,aes(x=date))+geom_line(aes(y=E,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Exposed (count)", labels = scales::comma)+labs(x="Date",size=1.3)+
  geom_line(mapping=aes(y=Ehat*N,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank()) 
p3=ggplot(df,aes(x=date))+geom_line(aes(y=I,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Infectious (count)", labels = scales::comma)+labs(x="Date",size=1.3)+
  geom_line(mapping=aes(y=Ihat*N,col="Fitted"),size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank()) 
p4=ggplot(df,aes(x=date))+geom_line(aes(y=R,col="Observed"),size=1.3)+
  theme(axis.text.x = element_text(angle = 90))+scale_x_date(date_labels = "%Y.%b.%d",limits =c(start,end))+
  geom_vline(xintercept = as_date(c('2021-04-08')), color = 'orange', size = 1.3,linetype="dotted") +
  geom_vline(xintercept = as_date(c('2020-12-26')), color = 'orange', size = 1.3,linetype="dotted") +
  scale_y_continuous(name="Recovered (count)", labels = scales::comma)+
  geom_line(mapping=aes(y=Rhat*N,col="Fitted"),size=1.3)+labs(x="Date",size=1.3)+
  theme(axis.text=element_text(size=7))+ theme(legend.title = element_blank()) 

ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

cases_v2=read.csv(file.path(directory, "data/covid19-download_v2.csv"))
#cases <- read.csv("covid19-download.csv",header=T)
mydata2 = cases_v2 %>% filter(pruid==35) %>%
  select(prname,date,numtoday,numdeathstoday,numrecoveredtoday,numactive,numtestedtoday,numteststoday,numtests,numrecover) %>%
  rowwise() %>% mutate(numtestsmax=max(numtestedtoday,numteststoday,na.rm=T))%>%
  select(-numtestedtoday,-numteststoday) %>%
  drop_na() %>% mutate(date=as_date(date)) %>%
  mutate(conftotal=cumsum(numtoday)) %>%
  mutate(N=14750653) %>%
  mutate(E=rollapplyr(numtoday,width=latent_period,FUN=sum,partial=T)) %>%
  mutate(I=numactive-E ,R=numrecover) %>% drop_na(conftotal)%>%
  complete(date = full_seq(date, period = 1)) %>% mutate(S=N-E-I-R) 
start=mydata2$date[1]
mydata2$time=mydata2$date-start+1
df2=mydata2[which(mydata2$date==as_date("2021-05-21")):length(mydata2$date),]
df2=df2[1:interval,]
beta=mydata[which(mydata$date==as_date("2021-05-21")),c("beta")][[1]] %>%as.numeric
output_v2=seir(beta,gamma,delta,df2$S[1],df2$E[1],df2$I[1],df2$R[1],1) %>%as.data.frame()

#
par(mfrow=c(2, 2))
N=14750653
ymax=max(df2$S,output_v2$S*N,na.rm=T)
ymin=min(df2$S,output_v2$S*N,na.rm=T)
#xmax=mydata$date[1]
#xmin=mydata$date[length(mydata$date)]
matplot(df2$time,df2$S,type = 'l',col='blue',ylim=c(ymin,ymax),xlab="time",ylab="count",main="Susceptible")
par(new = TRUE)
matplot(output_v2$time,output_v2$S*N,type = 'l',col='red',ylim=c(ymin,ymax),xlab="time",ylab="count")

ymax=max(df2$E,output_v2$E*N,na.rm=T)
ymin=min(df2$E,output_v2$E*N,na.rm=T)
matplot(df2$time,df2$E,type = 'l',col='blue',ylim=c(ymin,ymax),xlab="time",ylab="count",main="Exposed")
par(new = TRUE)
matplot(output_v2$time,output_v2$E*N,type = 'l',col='red',ylim=c(ymin,ymax),xlab="time",ylab="count")

ymax=max(df2$I,output_v2$I*N,na.rm=T)
ymin=min(df2$I,output_v2$I*N,na.rm=T)
matplot(df2$time,df2$I,type = 'l',col='blue',ylim=c(ymin,ymax),xlab="time",ylab="count",main="Infectious")
par(new = TRUE)
matplot(output_v2$time,output_v2$I*N,type = 'l',col='red',ylim=c(ymin,ymax),xlab="time",ylab="count")

ymax=max(df2$R,output_v2$R*N,na.rm=T)
ymin=min(df2$R,output_v2$R*N,na.rm=T)
matplot(df2$time,df2$R,type = 'l',col='blue',ylim=c(ymin,ymax),xlab="time",ylab="count",main="Recovered")
par(new = TRUE)
matplot(output_v2$time,output_v2$R*N,type = 'l',col='red',ylim=c(ymin,ymax),xlab="time",ylab="count")


## Codes updated by Hana F.
png('Predictions.png', width = 600, height = 500, pointsize = 19)
par(mfrow=c(2, 2), oma=c(2,0,0,0), mar=c(3,3.5,2,1), mgp=c(2.75,0.8,0), font.lab=2, 
    cex.axis=0.8, cex.lab=0.8, cex.main=1.0)
N=14750653

byMillions <- 1000000 # for scaling y-axis
byThousands <- 1000 # for scaling y-axis

# Combine df2 with output_v2 to df3
names(output_v2) <- paste0(names(output_v2), '_est')
df3 <- cbind(df2, output_v2)
df3$S_N <- df3$S_est*N
df3$E_N <- df3$E_est*N
df3$I_N <- df3$I_est*N
df3$R_N <- df3$R_est*N

# ymax=max(df2$S,output_v2$S*N,na.rm=T)/byMillions
# ymin=min(df2$S,output_v2$S*N,na.rm=T)/byMillions
#xmax=mydata$date[1]
#xmin=mydata$date[length(mydata$date)]
ymax=max(df3[c('S','S_N')],na.rm=T)/byMillions
ymin=min(df3[c('S','S_N')],na.rm=T)/byMillions

dateNum <- as.numeric(df3$date) # for setting x-axis labeling
xAt <- seq(dateNum[1], tail(dateNum,1), by=5) # for setting x-axis labeling

matplot(df3$date,df3[c('S','S_N')]/byMillions,type = 'l', lty = 1, lwd = 1.5,
        xaxt = 'n',yaxt = 'n',
        col=c('blue','red'), ylim=c(ymin,ymax),
        xlab="", ylab="Population (millions)", main="Susceptible")
text(x = mean(dateNum), y = ymin - 0.02, labels = "Date", xpd = NA, srt = 0, cex = 0.85, font = 2)
axis(1, at=xAt, labels = FALSE)
text(x= xAt, y = par("usr")[3] - 0.008,
     labels = format(as.Date(xAt), '%B %d'), 
     xpd = NA,
     srt = 45, ## Rotate the labels by 45 degrees.
     cex = 0.75)
axis(2, las = 2)
# par(new = TRUE)
# matplot(output_v2$time,output_v2$S*N,type = 'l',col='red',ylim=c(ymin,ymax),xlab="time",ylab="count")

ymax=max(df3[c('E','E_N')],na.rm=T)
ymin=min(df3[c('E','E_N')],na.rm=T)

matplot(df3$date,df3[c('E','E_N')],type = 'l', lty = 1, lwd = 1.5,
        xaxt = 'n',yaxt = 'n',
        col=c('blue','red'), ylim=c(ymin,ymax),
        xlab="", ylab="Population", main="Exposed")
text(x = mean(dateNum), y = ymin - 3500, labels = "Date", xpd = NA, srt = 0, cex = 0.85, font = 2)
axis(1, at=xAt, labels = FALSE)
text(x= xAt, y = par("usr")[3] - 1500,
     labels = format(as.Date(xAt), '%B %d'), 
     xpd = NA,
     srt = 45, ## Rotate the labels by 45 degrees.
     cex = 0.75)
axis(2, las = 2)

ymax=max(df3[c('I','I_N')],na.rm=T)
ymin=min(df3[c('I','I_N')],na.rm=T)

matplot(df3$date,df3[c('I','I_N')],type = 'l', lty = 1, lwd = 1.5,
        xaxt = 'n',yaxt = 'n',
        col=c('blue','red'), ylim=c(ymin,ymax),
        xlab="", ylab="Population", main="Infectious")
text(x = mean(dateNum), y = ymin - 5500, labels = "Date", xpd = NA, srt = 0, cex = 0.85, font = 2)
axis(1, at=xAt, labels = FALSE)
text(x= xAt, y = par("usr")[3] - 2300,
     labels = format(as.Date(xAt), '%B %d'), 
     xpd = NA,
     srt = 45, ## Rotate the labels by 45 degrees.
     cex = 0.75)
axis(2, las = 2)

ymax=max(df3[c('R','R_N')],na.rm=T)/byThousands
ymin=min(df3[c('R','R_N')],na.rm=T)/byThousands

matplot(df3$date,df3[c('R','R_N')]/byThousands,type = 'l', lty = 1, lwd = 1.5,
        xaxt = 'n',yaxt = 'n',
        col=c('blue','red'), ylim=c(ymin,ymax),
        xlab="", ylab="Population (thousands)", main="Recovered")
text(x = mean(dateNum), y = ymin - 20, labels = "Date", xpd = NA, srt = 0, cex = 0.85, font = 2)
axis(1, at=xAt, labels = FALSE)
text(x= xAt, y = par("usr")[3] -8.5,
     labels = format(as.Date(xAt), '%B %d'), 
     xpd = NA,
     srt = 45, ## Rotate the labels by 45 degrees.
     cex = 0.75)
axis(2, las = 2)

# Legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Observed","Estimated"), xpd = TRUE, horiz = TRUE, 
       inset = c(0,0), bty = "n", lty = 1, lwd=2, col = c("blue","red"), cex = 0.9)

dev.off()


# RMSE calculation

observed <- df2[, c('time', 'S', 'E', 'I', 'R')]
predicted <- output_v2 %>% mutate(S = S*N, E = E*N, I = I*N, R=R*N)

rmse.s <- rmse(observed$S, predicted$S) #818677.3
rmse.e <- rmse(observed$E, predicted$E) #3651.17
rmse.i <- rmse(observed$I, predicted$I) #7948.69
rmse.r <- rmse(observed$R, predicted$R) #30716.83

#relative RMSE

rrmse.s <- rmse.s/sd(observed$S) #123.8984
rrmse.e <- rmse.e/sd(observed$E) #1.480636
rrmse.i <- rmse.i/sd(observed$I) #2.323698
rrmse.r <- rmse.r/sd(observed$R) #2.466726
