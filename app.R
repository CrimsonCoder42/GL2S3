#GL2Z3
library(shiny)
library(shinydashboard)
library(shinyWidgets)
source("jaxmat.R")
#The user interface
header <- dashboardHeader(
  title = HTML("GL(2,Z<sub>3</sub>)"),
  titleWidth = 600
)
sidebar <- dashboardSidebar(disable = TRUE)
body <- dashboardBody(
  fluidRow(
    column(
      width = 3,
      title = "Create a Matrix",
      radioButtons("det", "Choose a determinant",
                   choiceNames = c("Determinant 1","Determinant -1"),
                   choiceValues = c(1,-1)
      ),
      radioButtons("trace", "Choose a trace",
        choiceNames = c("Trace 0","Trace 1","Trace -1"),
        choiceValues = c(0,1,-1)
      ),
      actionButton("generate","Create a Matrix"),
      uiOutput("matrixA"),
      actionButton("applysub","Apply Subspaces"),
      # actionButton("permute","Construct the Permutation"),
      actionButton("powers","Calculate Powers"),
      actionButton("multm1","Multiply A by -1"),
      # actionButton("mult2","Multiply A by 2")
    ),
    column(
      width = 3,
      h3("The Six Subspaces"),
      h4("Subspace 1"),
      h5(jax.vecList(rbind(c(1,1,2),c(0,0,0)))),
      h4("Subspace 2"),
      h5(jax.vecList(rbind(c(0,0,0),c(0,1,2)))),
      h4("Subspace 3"),
      h5(jax.vecList(rbind(c(0,1,2),c(0,1,2)))),
      h4("Subspace 4"),
      h5(jax.vecList(rbind(c(0,1,2),c(0,2,1))))
    ),
    column(
      width = 3,
      h3("Action on Subspaces"),
      h4("Action on Subspace 1"),
      uiOutput("sub1"),
      h4("Action on Subspace 2"),
      uiOutput("sub2"),
      h4("Action on Subspace 3"),
      uiOutput("sub3"),
      h4("Action on Subspace 4"),
      uiOutput("sub4"),

      h3("Cycle Representation"),
      uiOutput("perm")
         
    ),
    column(
      width = 3,
      h3("The Powers of A"),
      uiOutput("power1"),
      uiOutput("power2"),
      uiOutput("power3"),
      uiOutput("power4"),
      uiOutput("power5"),
      uiOutput("power6"),
      uiOutput("power7"),
      uiOutput("power8"),
      uiOutput("power9"),
      uiOutput("power10"),
      uiOutput("power11"),
      uiOutput("power12"),
      uiOutput("power16"),
    )
  )
)
ui <- dashboardPage(header, sidebar, body)

#Functions that implement the mathematics
source("Z5calc.R")
source("permutecalc.R")

#Functions that read the input and modify the output and input
server <- function(session, input, output) {
  #Variables that are shared among server functions
  A <- matrix(nrow = 2, ncol = 2)
  
  #Functions that respond to events in the input
  #According to the documentation, a global withMathJax() will not work here
  observeEvent(input$generate,{
     A <<- Z3CreateMatrix(as.numeric(input$det),as.numeric(input$trace))
     output$matrixA <- renderUI({jax.matrix(A, name = "A")})
     output$sub1 <- renderUI("")
     output$sub2 <- renderUI("")
     output$sub3 <- renderUI("")
     output$sub4 <- renderUI("")
     output$sub5 <- renderUI("")
     output$sub6 <- renderUI("")
     output$power1 <- renderUI("")
     output$power2 <- renderUI("")
     output$power3 <- renderUI("")
     output$power4 <- renderUI("")
     output$power5 <- renderUI("")
     output$power6 <- renderUI("")
     output$perm <- renderUI("")
  })
  observeEvent(input$applysub,{
      v1 <- c(1,0)
      x <- ActOnVectorinZ3(A, v1)
      output$sub1 <- renderUI({jax.mTimesV("A",v1,x)})
      v2 <- c(0,1)
      x2 <- ActOnVectorinZ3(A, v2)
      output$sub2 <- renderUI({jax.mTimesV("A",v2,x2)})
      v3 <- c(1,1)
      x3 <- ActOnVectorinZ3(A, v3)
      output$sub3 <- renderUI({jax.mTimesV("A",v3,x3)})
      v4 <- c(1,-1)
      x4 <- ActOnVectorinZ3(A, v4)
      output$sub4 <- renderUI({jax.mTimesV("A",v4,x4)})
  })
  observeEvent(input$permute,{
      fval <- sapply(1:6,TransforminZ3,A=A)
      output$perm <-renderUI(h3(Perm.cycle.convert(fval)))
    
  })
#It may be possible to put the MathJax stuff into a function
  observeEvent(input$powers,{
      output$power1 <- renderUI({jax.matrix(A, name = "A")})
      A2 <- Z3MatProd(A,A)
      output$power2 <- renderUI({jax.matrix(A2, name = "A^2")})
      
      A3 <- Z3MatProd(A,A2)
      output$power3 <- renderUI({jax.matrix(A3, name = "A^3")})
      
      A4 <- Z3MatProd(A,A3)
      output$power4 <- renderUI({jax.matrix(A4, name = "A^4")})
      
      A5 <- Z3MatProd(A,A4)
      output$power5 <- renderUI({jax.matrix(A5, name = "A^5")})
      
      A6 <- Z3MatProd(A,A5)
      output$power6 <- renderUI({jax.matrix(A6, name = "A^6")})
      
      A7 <- Z3MatProd(A,A6)
      output$power7 <- renderUI({jax.matrix(A7, name = "A^7")})
      
      A8 <- Z3MatProd(A,A7)
      output$power8 <- renderUI({jax.matrix(A8, name = "A^8")})
      
      A9 <- Z3MatProd(A,A8)
      output$power9 <- renderUI({jax.matrix(A9, name = "A^9")})
      
      A10 <- Z3MatProd(A,A9)
      output$power10 <- renderUI({jax.matrix(A10, name = "A^{10}")})
      
      A11 <- Z3MatProd(A,A10)
      output$power11 <- renderUI({jax.matrix(A11, name = "A^{11}")})
      
      A12 <- Z3MatProd(A,A11)
      output$power12 <- renderUI({jax.matrix(A12, name = "A^{12}")})
      
      A13 <- Z3MatProd(A,A12)
      A14 <- Z3MatProd(A,A13)
      A15 <- Z3MatProd(A,A14)
      A16 <- Z3MatProd(A,A15)
      output$power16 <- renderUI({jax.matrix(A16, name = "A^{16}")})
  })

  #Multiply by -1
  observeEvent(input$multm1,{
    A <<- matrix(vZ5Prod(A,-1),nrow = 2)
    output$matrixA <- renderUI({jax.matrix(A, name = "A")})
  })

  #Multiply by 2
observeEvent(input$mult2,{
  A <<- matrix(vZ5Prod(A,2),nrow = 2)
  output$matrixA <- renderUI({jax.matrix(A, name = "A")})
  })
}
  

#Run the app
shinyApp(ui = ui, server = server)