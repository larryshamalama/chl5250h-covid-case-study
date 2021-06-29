library(shiny)
library(dplyr)
library(lubridate)
library(readr)
library(deSolve)
library(rstudioapi)
library(zoo)
library(tidyr)
library(Metrics) 

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

#function to calculated different measures of error
rmse <- function(pred, obs){
  mse <- mean((pred - obs)^2, na.rm = TRUE)
  rmse <- sqrt(mse)
  rrmse <- rmse/sd(obs, na.rm = TRUE)
  c(mse, rmse, rrmse)
}

# Parameters
t=7

latent_period <- 5
interval=30
infectious_period <- 9

gamma <- 1/infectious_period
delta <- 1/latent_period

N=14750653

# Load RDS
df2 <- readRDS('data/df2.rds')
observed <- df2[, c('time', 'S', 'E', 'I', 'R')]

# ShinyApp
ui <- fluidPage(
   
   
   titlePanel("Find Best Match Beta"),
     
      # mainPanel(
      #   sliderInput("beta",
      #               "Value of Beta:",
      #               min = 0.01,
      #               max = 0.4,
      #               value = 0.1145993),
      #   plotOutput("EsPlot",  
      #              width = 700,
      #              height = 600),
      #   tableOutput("rmseTable")
      # )
   
   fluidRow(
     column(width = 6.5,
            sliderInput("beta",
                        "Value of Beta:",
                        min = 0.01,
                        max = 0.4,
                        value = 0.1145993),
            plotOutput("EsPlot",
                       width = 700,
                       height = 600)
    ),
    column(width = 4,
           tableOutput("rmseTable"))
   )
)

server <- function(input, output) {
  
  output_v2 <- reactive({
    seir(input$beta,gamma,delta,df2$S[1],df2$E[1],df2$I[1],df2$R[1],1) %>% as.data.frame()
  })
  
  predicted <- reactive({
    output_v2() %>% mutate(S = S*N, E = E*N, I = I*N, R=R*N)
  })
  
  output$EsPlot <- renderPlot({
    par(mfrow=c(2, 2), oma=c(4,0,0,0), mar=c(4,4.5,2,1), mgp=c(3.5,0.8,0), cex = 1.2, font.lab=2)
    
    byMillions <- 1000000 # for scaling y-axis
    byThousands <- 1000 # for scaling y-axis
    
    # Combine df2 with output_v2 to df3
    output_v2b <- output_v2()
    names(output_v2b) <- paste0(names(output_v2b), '_est')
    df3 <- cbind(df2, output_v2b)
    df3$S_N <- df3$S_est*N
    df3$E_N <- df3$E_est*N
    df3$I_N <- df3$I_est*N
    df3$R_N <- df3$R_est*N
    
    ymax=max(df3[c('S','S_N')],na.rm=T)/byMillions
    ymin=min(df3[c('S','S_N')],na.rm=T)/byMillions
    
    dateNum <- as.numeric(df3$date) # for setting x-axis labeling
    xAt <- seq(dateNum[1], tail(dateNum,1), by=5) # for setting x-axis labeling
    
    matplot(df3$date,df3[c('S','S_N')]/byMillions,type = 'l', lty = 1, lwd = 1.5,
            xaxt = 'n',yaxt = 'n',
            col=c('blue','red'), ylim=c(ymin,ymax),
            xlab="", ylab="Population (millions)", main="Susceptible")
    text(x = mean(dateNum), y = ymin - (ymax-ymin)/2.5, labels = "Date", xpd = NA, srt = 0, cex = 1.1, font = 2)
    axis(1, at=xAt, labels = FALSE)
    text(x= xAt, y = par("usr")[3] - (ymax-ymin)/6,
         labels = format(as.Date(xAt), '%B %d'), 
         xpd = NA,
         srt = 45, ## Rotate the labels by 45 degrees.
         cex = 0.9)
    axis(2, las = 2)

    ymax=max(df3[c('E','E_N')],na.rm=T)
    ymin=min(df3[c('E','E_N')],na.rm=T)
    
    matplot(df3$date,df3[c('E','E_N')],type = 'l', lty = 1, lwd = 1.5,
            xaxt = 'n',yaxt = 'n',
            col=c('blue','red'), ylim=c(ymin,ymax),
            xlab="", ylab="Population", main="Exposed")
    text(x = mean(dateNum), y = ymin - (ymax-ymin)/2.5, labels = "Date", xpd = NA, srt = 0, cex = 1.1, font = 2)
    axis(1, at=xAt, labels = FALSE)
    text(x= xAt, y = par("usr")[3] - (ymax-ymin)/6,
         labels = format(as.Date(xAt), '%B %d'), 
         xpd = NA,
         srt = 45, ## Rotate the labels by 45 degrees.
         cex = 0.9)
    axis(2, las = 2)
    
    ymax=max(df3[c('I','I_N')],na.rm=T)
    ymin=min(df3[c('I','I_N')],na.rm=T)
    
    matplot(df3$date,df3[c('I','I_N')],type = 'l', lty = 1, lwd = 1.5,
            xaxt = 'n',yaxt = 'n',
            col=c('blue','red'), ylim=c(ymin,ymax),
            xlab="", ylab="Population", main="Infectious")
    text(x = mean(dateNum), y = ymin - (ymax-ymin)/2.5, labels = "Date", xpd = NA, srt = 0, cex = 1.1, font = 2)
    axis(1, at=xAt, labels = FALSE)
    text(x= xAt, y = par("usr")[3] - (ymax-ymin)/6,
         labels = format(as.Date(xAt), '%B %d'), 
         xpd = NA,
         srt = 45, ## Rotate the labels by 45 degrees.
         cex = 0.9)
    axis(2, las = 2)
    
    ymax=max(df3[c('R','R_N')],na.rm=T)/byThousands
    ymin=min(df3[c('R','R_N')],na.rm=T)/byThousands
    
    matplot(df3$date,df3[c('R','R_N')]/byThousands,type = 'l', lty = 1, lwd = 1.5,
            xaxt = 'n',yaxt = 'n',
            col=c('blue','red'), ylim=c(ymin,ymax),
            xlab="", ylab="Population (thousands)", main="Recovered")
    text(x = mean(dateNum), y = ymin - (ymax-ymin)/2.5, labels = "Date", xpd = NA, srt = 0, cex = 1.1, font = 2)
    axis(1, at=xAt, labels = FALSE)
    text(x= xAt, y = par("usr")[3] -(ymax-ymin)/6,
         labels = format(as.Date(xAt), '%B %d'), 
         xpd = NA,
         srt = 45, ## Rotate the labels by 45 degrees.
         cex = 0.9)
    axis(2, las = 2)
    
    # Legend
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", c("Observed","Estimated"), xpd = TRUE, horiz = TRUE, 
           inset = c(0,0.05), bty = "n", lty = 1, lwd=2, col = c("blue","red"), cex = 1.2)
  })
  

  output$rmseTable <- renderTable({
    rmse.s <- rmse(observed$S, predicted()$S)
    rmse.e <- rmse(observed$E, predicted()$E)
    rmse.i <- rmse(observed$I, predicted()$I)
    rmse.r <- rmse(observed$R, predicted()$R)
    
    rmseTable <- data.frame(SEIR=character(),
                            mse=double(),
                            rmse=double(),
                            rrmse=double())
    
    rmseTable[1,] <- c('Susceptible', rmse.s)
    rmseTable[2,] <- c('Exposed', rmse.e)
    rmseTable[3,] <- c('Infectious', rmse.i)
    rmseTable[4,] <- c('Recovered', rmse.r)
    
    return(rmseTable)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

