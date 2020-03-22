
library(shiny)
library(shinydashboard)
library(shinyLP)
library(seqinr)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
require(Biostrings)


ui <- fluidPage(

    navbarPage("GAtRS",
               tabPanel("Home",
                        jumbotron(paste("Hi Welcome to GAtRS"), paste("This is web application for Genome Analysis using R Shiny."), buttonLabel = "Click Me" ), 
                        tags$a(href="https://github.com/Kamakshaiah", "Know More here")
                        
                        ),
               tabPanel("Data",
                        sidebarLayout(
                            sidebarPanel(
                                fileInput("file1", "Choose Sequnce", accept=c('text/csv', 'text/comma-separated-values', 'text/plain', '.csv', '.fasta')),
                                fileInput("file2", "Choose Comparison Sequence", accept=c('text/csv', 'text/comma-separated-values', 'text/plain', '.csv', '.fasta')),
                                hr(),
                                helpText("Upload .fasta Files"),
                                hr(),
                                helpText("For sample data sets visit:"),
                                tags$a(href="https://www.ncbi.nlm.nih.gov/nuccore/MN988668", "COVID-19")
                        
                            ), 
                            mainPanel(
                                verbatimTextOutput('seq1'),
                                verbatimTextOutput("seq2")
                            )
                        )
               ),
               navbarMenu("Analysis", 
                          tabPanel("GC Content",
                                   sidebarLayout(
                                       sidebarPanel(
                                           radioButtons("method", "Select Method", c("GC Content", "Counts", "GC Plot", "Rho", "T Test", "Chi Square", "Mutul Information", "Association" ))
                                           
                                       ),
                                       mainPanel(
                                           div("Analysis",
                                               verbatimTextOutput("analout")
                                           ),
                                           div("Plot",
                                               plotOutput("seqplot")
                                           )
                                       )
                                   )
                          ),
                          tabPanel("Pairwise Sequence Alignment",
                                   sidebarLayout(
                                       sidebarPanel(
                                           radioButtons("pair", "Select Option", c("Info", "Dot Plot1", "PairWise Global alignments2", "Pairwise local alignment3", "Multinomial Model4")),
                                           tags$small(helpText("1 Not implemented")),
                                           tags$small(helpText("2 [Needleman-Wunsch algorithm] Only basic methods are allowed; For advanced upgrade to Advanced Application")),
                                           tags$small(helpText("3 [Smith-Waterman algorithm]; Upgrade to advanced application")),
                                           tags$small(helpText("4 Upgrade to advanced application"))
                                       ),
                                       mainPanel(
                                          div(
                                              verbatimTextOutput("pairanal")
                                          )
                                        )
                                   )
                                   ),
                          tabPanel("Multiple Alignment and Phylogenetic trees",
                                   sidebarLayout(
                                       sidebarPanel(
                                           radioButtons('multiopt', 'Select Option', c("Info", "Get Sequences", "Multiple Alignments", "Genetic Distances", "phylogenetic tree")),
                                           textInput('seqs', "Accessions", value = "", placeholder = "Write Accession Numbers sep. by Commas")
                                       ),
                                       mainPanel(
                                           div(
                                               textOutput("multianal")
                                           )
                                       )
                                   )

                                   ),
                          tabPanel("Computational Gene-finding",
                                   sidebarLayout(
                                       sidebarPanel(
                                           radioButtons("cgfopt", "Select Option", c("Get SGC 1", "Match Pattern", "Tranlation 2" )),
                                           textInput('pattern', "Pattern", value = "atg", placeholder = "Write Pattern"),
                                           hr(),
                                           tags$small(helpText("1 Standard Genetic Code")),
                                           tags$small(helpText("2 Not implemented. Upgrade to Advanced Application."))
                                       ),
                                       mainPanel(
                                           div(
                                               verbatimTextOutput("cgfanal")
                                           ),
                                           div(
                                               plotOutput("cgfplot")
                                           )
                                       )
                                   )
                                   ),
                          tabPanel("Comparative Genomics",
                                   sidebarLayout(
                                       sidebarPanel(
                                           radioButtons("compmethod", "Select Method", c("Info", "T Test", "Chi Square", "Mutul Information", "Association" ))
                                       ),
                                       mainPanel(
                                           textOutput("companal")
                                       )
                                   )
                                   ),
                          tabPanel("Hidden Markov Models",
                                   sidebarLayout(
                                       sidebarPanel(
                                           radioButtons('hmmopt', 'Select Method', c("Info", "Simulate DNA using HMM", "HMM emission matrix", "State Inference"))
                                       ),
                                       mainPanel(
                                           textOutput("hmmanal")
                                       )
                                   )
                                   )
                          
                          
                          ),
               tabPanel("Contact",
                        sidebarLayout(
                            sidebarPanel("Contact Information"),
                            mainPanel(
                                div(
                                    htmlOutput("contactInfo")
                                )
                            )
                        )
               )
               )
               
               


)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    data_input1 <- reactive({
        infile <- input$file1
        req(infile)
        read.fasta(infile$datapath)
        
    } 
    )
    
    data_input2 <- reactive({
        infile <- input$file2
        req(infile)
        read.fasta(infile$datapath)
    })
    
    
    output$seq1 <- renderPrint({
        df <- data_input1()
        print(df)
    })
    
    output$seq2 <- renderPrint({
        df <- data_input2()
        print(df)
    })
    
    slidingwindowplot <- function(windowsize, inputseq)
    {
        starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
        n <- length(starts)    # Find the length of the vector "starts"
        chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
        for (i in 1:n) {
            chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
            chunkGC <- GC(chunk)
            print(chunkGC)
            chunkGCs[i] <- chunkGC
        }
        plot(starts,chunkGCs,type="b", xlab="Nucleotide start position", ylab="GC content", main="Window Plot for GC content for 1000 segments")
    }
    
    methods <- reactive({
        
        dfone <- data_input1()
        dataone <- dfone[[1]]
        
        dftwo <- data_input2()
        datatwo <- dftwo[[1]]
        
        library(MASS)
        
        if (input$method == "GC Content"){
            return(GC(dataone))
        } else if (input$method == "Counts"){
            out <- list()
            for (i in 1:4){
                
                out[[i]] <- count(dataone, i)
            }
            return(out)
        } else if (input$method == "GC Plot"){
            slidingwindowplot(1000, dataone)
        } else if (input$method == "Rho"){
            countslet <- count(dataone, 1)
            countstwo <- count(dataone, 2)
            g <- countslet['g']/sum(countslet)
            c <- countslet['c']/sum(countslet)
            gc <- countstwo['gc']/sum(countstwo)
            return(gc/(g*c))
            
        } else if (input$method == "T Test"){
            
            
            tab1 <- table(dataone)
            tab2 <- table(datatwo)
            
           return(t.test(tab1, tab2)) # two sample equality test
            return(out)
        } else if (input$method == "Chi Square"){
            
            tab1 <- table(dataone)
            tab2 <- table(datatwo)
            
            out <- chisq.test(tab1, tab2) # two sample equality test
            return(out)
        } else if (input$method == "Mutul Information"){
            
            tab1 <- table(dataone)
            tab2 <- table(datatwo)
            
            mytable <- xtabs(~tab1+tab2, data=cbind(tab1, tab2)) 
            out <- loglm(~tab1+tab2, mytable)
            return(out)
        } else if (input$method == "Association"){
            tab1 <- table(dataone)
            tab2 <- table(datatwo)
            
            out <- list()
            out[[4]] <- tab1
            out[[5]] <- tab2 
            for (i in 1:3){
                methods = c("pearson","kendall","spearman")
                out[[i]] <- cor(tab1, tab2, method=methods[[i]])
                
            }
            return(out)
        }
        
        
    })
    
    output$analout <- renderPrint({
        print(methods())
    })
    
    output$seqplot <- renderPlot({
        if (input$method == "GC Plot"){
            df <- data_input1()
            data <- df[[1]]
            slidingwindowplot(1000, data)
        }
        
    })
    
    
    pairout <- reactive({
        
        if (input$pair == "Info"){
            text1 <- paste("Requires minimum 16 GB RAM. Upgrade to Advanced Application. Visit CONTACT menu.")
            return(text1)
        }
        
        
    #     dfone <- data_input1()
    #     dataone <- dfone[[1]]
    # 
    #     dftwo <- data_input2()
    #     datatwo <- dftwo[[1]]
    # 
    #     if (input$pair == "Dot Plot"){
    #         dotPlot(dataone, datatwo)
    #     }
    #     if (input$pair == "PairWise Global alignments2"){
    #         data(BLOSUM50)
    #         
    #         seq1 <- c2s(dataone)
    #         seq2 <- c2s(datatwo)
    #         seq1 <- toupper(seq1)
    #         seq2 <- toupper(seq2)
    #         
    #         out <- pairwiseAlignment(seq1, seq2,
    #                           substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
    #         return(out)
    #     }
    })
    
    output$pairanal <- renderText({
        print(pairout())
    })
    
    multiout <- reactive({
        if (input$multiopt == "Info"){
            text <- paste("Requires minimum 16 GB RAM. Upgrade to Advanced Application. Visit CONTACT menu.")
            return(text)
        }
    })
    
    output$multianal <- renderText({
        print(multiout())
    })
    
        
    cgfout <- reactive({
        dfone <- data_input1()
        dataone <- dfone[[1]]
        
        if (input$cgfopt == "Match Pattern"){
            ptrn <- input$pattern
            dataone <- c2s(dataone)
            return(matchPattern(ptrn, dataone))
        }
    })    
    
    output$cgfplot <- renderPlot({
        if(input$cgfopt == "Get SGC 1"){
         tablecode()   
        }
    })
    
    output$cgfanal <- renderPrint({
        print(cgfout())
    })
    
    compout <- reactive({
        if (input$compmethod == "Info"){
            text <- paste("Go to menu 'GC Content'.")
            print(text)
        }
    })
    
    output$companal <- renderText({
        compout()
    })
    
    hmmout <- reactive({
        if (input$hmmopt == "Info"){
            text <- paste("Upgrade to Advanced Application. Get information from CONTACT menu.")
            print(text)
        }
    })
    
    output$hmmanal <- renderText({
        hmmout()
    })
    
    
    output$contactInfo <- renderText({
        text1 <- paste("Dr. M. Kamakshaiah")
        text2 <- paste("+919848396972")
        text3 <- paste("kamakshaiah.m@gmail.com")
        HTML(paste(text1, text2, text3, sep = '<br/>'))
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
