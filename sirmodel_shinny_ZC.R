library(shiny)
library(dplyr)
library(lubridate)
library(readr)
library(deSolve)
library(rstudioapi)
library(zoo)
library(tidyr)
library(Metrics) 


ui <- fluidPage(
   
   
   titlePanel("Find Best Match Beta"),
     
      mainPanel(
        sliderInput("beta",
                    "Value of Beta:",
                    min = 0.01,
                    max = 0.4,
                    value = 0.1145993),
        plotOutput("EsPlot"),
        h4("relative RMSE of Susceptible"),
        textOutput("rrmse.s"),
        h4("relative RMSE of Exposed"),
        textOutput("rrmse.e"),
        h4("relative RMSE of Infected"),
        textOutput("rrmse.i"),
        h4("relative RMSE of Recovered"),
        textOutput("rrmse.r")
        
      )
   )



server <- function(input, output) {
  
  output$EsPlot <- renderPlot({
    load("data/sumin_image.RData")
    output_v2=seir(input$beta,gamma,delta,df2$S[1],df2$E[1],df2$I[1],df2$R[1],1) %>% as.data.frame()
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
    matplot(df2$time,df2$I,type = 'l',col='blue',ylim=c(ymin,ymax),xlab="time",ylab="count",main="Infected")
    par(new = TRUE)
    matplot(output_v2$time,output_v2$I*N,type = 'l',col='red',ylim=c(ymin,ymax),xlab="time",ylab="count")
    
    ymax=max(df2$R,output_v2$R*N,na.rm=T)
    ymin=min(df2$R,output_v2$R*N,na.rm=T)
    matplot(df2$time,df2$R,type = 'l',col='blue',ylim=c(ymin,ymax),xlab="time",ylab="count",main="Recovered")
    par(new = TRUE)
    matplot(output_v2$time,output_v2$R*N,type = 'l',col='red',ylim=c(ymin,ymax),xlab="time",ylab="count")
  })
  
    
  
  output$rrmse.s <- renderText({
    load("data/sumin_image.RData")
    output_v2=seir(input$beta,gamma,delta,df2$S[1],df2$E[1],df2$I[1],df2$R[1],1) %>%as.data.frame()
    predicted <- output_v2 %>% mutate(S = S*N, E = E*N, I = I*N, R=R*N)
    rmse.s <- rmse(observed$S, predicted$S)
    rmse.s/sd(observed$S)
  })
  
  output$rrmse.e <- renderText({
    load("data/sumin_image.RData")
    output_v2=seir(input$beta,gamma,delta,df2$S[1],df2$E[1],df2$I[1],df2$R[1],1) %>%as.data.frame()
    predicted <- output_v2 %>% mutate(S = S*N, E = E*N, I = I*N, R=R*N)
    rmse.e <- rmse(observed$E, predicted$E)
    rmse.e/sd(observed$E)
  })
  
  output$rrmse.i <- renderText({
    load("data/sumin_image.RData")
    output_v2=seir(input$beta,gamma,delta,df2$S[1],df2$E[1],df2$I[1],df2$R[1],1) %>% as.data.frame()
    predicted <- output_v2 %>% mutate(S = S*N, E = E*N, I = I*N, R=R*N)
    rmse.i <- rmse(observed$I, predicted$I)
    rmse.i/sd(observed$I)
  })
  
  output$rrmse.r <- renderText({
    load("data/sumin_image.RData")
    output_v2=seir(input$beta,gamma,delta,df2$S[1],df2$E[1],df2$I[1],df2$R[1],1) %>% as.data.frame()
    predicted <- output_v2 %>% mutate(S = S*N, E = E*N, I = I*N, R=R*N)
    rmse.r <- rmse(observed$R, predicted$R)
    rmse.r/sd(observed$R)
  })
}


# Run the application 
shinyApp(ui = ui, server = server)

