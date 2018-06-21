library(shiny)
library(shinysky)
library(DT)

options(shiny.maxRequestSize=1000*1024^2)

ui <- fluidPage(
  titlePanel("AA Finder"),
  sidebarLayout(
    sidebarPanel(
      wellPanel(selectInput(inputId = "protein", label = "Choose a protein file",c("Protein 2015 / RELEASE.106" = "2015",
                                                "Protein 2016 / RELEASE.107" = "2016",
                                                "Protein 2018 / RELEASE.108" = "2018")),
      fileInput(inputId = "altProtein", label = "Alternatively, upload a protein file (.fa extention)",accept=c('.fasta','.fa'))),
      fileInput(inputId = "inputCsv", label = "Input CSV file",accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
      actionButton("sample", "Sample Data"),
      actionButton("submit", "Start Processing"),
      downloadButton("download", "Download Results"),
      busyIndicator("Preparing data",wait = 0),
      textOutput("done")),
  mainPanel(
    fluidRow(tags$p(strong("AA Finder is a tool which calculates the position of certain amino acids (Serine, Threonine, Tyrosine) within a larger protein based on its position in a peptide.")),
            tags$i("Usage:"), br(),
            "1. Choose the protein database to calculate from (All databases taken from ftp.ncbi.nih.gov) or upload a .fa file database", br(),
            "2. Input a CSV input file (Requires certain heading names: Gene, GI_accession, Peptide, Residue.Both)", br(),
            "3. Click 'Start Processing' ", br(), br(),
            tags$p("Any questions can be directed to",tags$a("Zilong Huang", href = "mailto:zilonghuang6@gmail.com?subject=AA%20Finder")," or ",tags$a("Gurkan Bebek", href = "mailto:gurkan@case.edu?subject=AA%20Finder"))),br(),
    textOutput("unmatched"),br(),
    fluidRow(DT::dataTableOutput("table"))))
)
AAFinder <- function(csv,protein,altprotein)
{
  library(seqinr)
  Input <- read.csv(csv, stringsAsFactors = FALSE)  #Input file must be named Input.csv
  data <- read.fasta(paste("protein",protein,".fa",sep = ""), as.string = TRUE)  #Protein database must be named protein.fa
  if(is.null(altprotein) == FALSE)
    data <- read.fasta(paste(altprotein,sep = ""), as.string = TRUE)
  
  if(length(grep("GI", names(Input))) != 0)
  {
    GI<- Input[ ,grep("GI", names(Input))]
  } else {
      GI<- Input[ ,3]
  }
  GIOrig <- GI
  GI <- gsub("-[0-9]","", GI)
  
  if(length(grep("Residue.Both",names(Input))) != 0)
  {
    AA <- Input["Residue.Both"][ ,1]
  } else if(length(grep("AA",names(Input))) != 0)
  {
    AA <- Input[ ,grep("AA",names(Input))]
  } else {
    AA <- Input[ , 5]
  }
  
  Unmatched <- 0
  Length <- length(GI)  #Finds how many values being searched for
  Peptide <- Input["Peptide"][ ,1]
  Peptide <- gsub('[a-z]',"",Peptide) #There were values such as "Mo" that would only show up on the protein as "M"
  
  StartingLetter <- gsub("\\d", "", AA) #Takes the letters
  StartingLetter <- gsub(";", "", StartingLetter) #These are for managing semicolons, extra commas, spaces and other extras
  StartingLetter <- gsub("  ", " ", StartingLetter)
  StartingLetter <- gsub(" ", ",", StartingLetter)
  StartingLetter <- gsub(",$", "", StartingLetter)
  StartingLetter <- gsub(", ", " ", StartingLetter)
  StartingLetter <- gsub("^,", " ", StartingLetter)
  LocationOnPeptide <- gsub("[A-z]", "", AA) #Takes the numbers
  LocationOnPeptide <- gsub(";", "", LocationOnPeptide) #Also managing extras
  LocationOnPeptide <- gsub("  ", " ", LocationOnPeptide)
  LocationOnPeptide <- gsub(" ", ",", LocationOnPeptide)
  LocationOnPeptide <- gsub(",$", "", LocationOnPeptide)
  LocationOnPeptide <- gsub(", ", " ", LocationOnPeptide)
  LocationOnPeptide <- gsub("^,", " ", LocationOnPeptide)
  
  OutputMatrix <- matrix(data = "", nrow=Length, ncol = 5)
  
  OutputMatrix[ , 1] <- Input["Gene"][ ,1] #Transfers values that stay constant
  OutputMatrix[ , 2] <- GIOrig
  OutputMatrix[ , 3] <- Input["Peptide"][ ,1]
  OutputMatrix[ , 4] <- AA
  
  SpacedOutMatrix <- matrix(data = " ", nrow = 6, ncol = Length)
  SpacedOutMatrix2 <- matrix(data = " ", nrow = 6, ncol = Length)
  S <- "S"  #These only exist because in SpacedOut2, the letters register as objects and not as characters
  T <- "T"
  Y <- "Y"
  dim(LocationOnPeptide) <- c(1,Length)
  dim(StartingLetter) <- c(1,Length)
  count3 <- 1
  while(count3 <= Length) { #This while loop spaces out values separated by commas
    if(grepl(",",LocationOnPeptide[1,count3]) == TRUE) {
      SpacedOut <- eval(parse(text = paste("list(", LocationOnPeptide[1,count3], ")"))) #Separates the comma'd values into a list
      SpacedOutMatrix[1:length(SpacedOut) ,count3] <- SpacedOut #Creates a matrix for the separated values
      dim(SpacedOutMatrix) <- c(6,Length)
      SpacedOut2 <- eval(parse(text=paste('list(', StartingLetter[1,count3], ')')))
      SpacedOutMatrix2[1:length(SpacedOut2) ,count3] <- SpacedOut2
      dim(SpacedOutMatrix2) <- c(6,Length)
      count3 <- count3 + 1
    }
    else  { count3 <- count3 + 1
    }
  }
  LocationOnProtein <- matrix(data = 1:Length, nrow = 1, ncol = Length)
  count <- 1 #Repeats finding the protein and pasting the exact protein in the ProteinMatrix
  ProteinMatrix <- matrix(data = 1:Length, nrow = 1, ncol= Length)
  while(count <= Length) {
    ProteinLocation <- grep(paste("|",GI[count],"|", sep = ""),names(data), fixed = TRUE)
    if(length(ProteinLocation) < 1 ) {
      ProteinMatrix[1,count] <- NA  #So if the grep value is integer(0), the length is 0 and put down as NA
      count <- count + 1
    }
    else {
      ProteinMatrix[1,count] <- ProteinLocation #Otherwise, values are just put into the matrix
      count <- count + 1
    }
  }
  LocationList <- matrix(data = "", nrow = 6, ncol = Length)
  dim(LocationList) <- c(6, Length)
  dim(SpacedOutMatrix2) <- c(6, Length)
  count2 <- 1
  while(count2 <= Length) {
    if(grepl(",",LocationOnPeptide[1,count2]) == TRUE) {  #If there is a comma on the location value
      count2 <- count2 + 1
      count4 <- 1
      while(count4 <= 6 && SpacedOutMatrix[count4,count2 - 1] != " ") {  #This records the location of every value
        LocationList[count4,count2 -1]  <- regexpr(Peptide[count2 -1], toupper(data[as.integer(ProteinMatrix[1,count2 -1])])) + as.integer(paste(unlist(SpacedOutMatrix[count4,count2 -1])), collapse='') -1
        if(regexpr(Peptide[count2 -1], toupper(data[as.integer(ProteinMatrix[1,count2 -1])])) <= 0)  {
          count4 <- 7
          OutputMatrix[count2-1, 5] <- "NoMatchFound" #If the peptide doesn't match, than there is no match found. If the GI didn't match, then the NA value still gives -1
          Unmatched <- Unmatched + 1
        } else {
          count4 <- count4 + 1  #If the match is found, it gets placed into the output matrix
          OutputMatrix[count2 -1, 5] <- paste(unlist(paste(SpacedOutMatrix2[,count2 -1],LocationList[ , count2 -1], sep = "")), collapse=",")
        }
      }
    }
    else  {
      LocationOnProtein[1, count2]  <- regexpr(Peptide[count2], toupper(data[as.integer(ProteinMatrix[1,count2])])) + as.integer(LocationOnPeptide[1,count2]) -1
      if(regexpr(Peptide[count2], toupper(data[as.integer(ProteinMatrix[1,count2])])) <= 0)  {
        OutputMatrix[count2, 5] <- "NoMatchFound" #If the peptide doesn't match
        Unmatched <- Unmatched + 1
      } else  {
        OutputMatrix[count2, 5] <- (paste(StartingLetter[1,count2],LocationOnProtein[1,count2], sep = ""))
      }
      count2 <- count2 + 1
    }
  }
  OutputMatrix[ ,5] <- gsub(" ," ," ", OutputMatrix[ ,5])
  OutputMatrix[ ,5] <- gsub(", " ," ", OutputMatrix[ ,5])
  UnmatchedNum <<- Unmatched
  Unmatched <<- paste("There were a total of ",Unmatched," Peptide/GI combinations that could not be found in the database. There were ", Length-Unmatched, " successfully processed entries")
  colnames(OutputMatrix) <- c("Gene","GI","Peptide","Peptide Location","Protein Location")
  return(OutputMatrix)
}

server <- function(input, output) {

  observeEvent(input$sample,
      {
        output$table <- DT::renderDataTable(DT::datatable({
          result <<-AAFinder("sampleInput.csv","2015",NULL)
        }))
        output$done <- renderText(paste("Processing completed at",Sys.time(),"using 2015 database with sample input"))
        output$unmatched <- renderText(Unmatched)
      
      output$download <- downloadHandler(
        filename = function(){
          paste("AAFinder_SAMPLE_(Unmatched: ",UnmatchedNum,")_",Sys.time(),".csv",sep = "")
        },
        content = function(file){
          write.csv(result, file)
        }
      )
      })
  
  
  observeEvent(input$submit, {
      if(is.null(input$inputCsv))
      {
        return()
      }
    
      output$table <- DT::renderDataTable(DT::datatable({
      req(input$inputCsv)
      result <<-AAFinder(input$inputCsv$datapath,input$protein,input$altProtein$datapath)
      }))
      
      output$unmatched <- renderText(Unmatched)
      
      output$download <- downloadHandler(
        filename = function(){
          paste("AAFinder_processed_(Unmatched: ",UnmatchedNum,")_",Sys.time(),".csv",sep = "")
        },
        content = function(file){
          write.csv(result, file)
        }
      )
      if(is.null(input$altProtein) == FALSE)
        output$done <- renderText(paste("Processing completed at",Sys.time(),"using user submitted data"))
      else
        output$done <- renderText(paste("Processing completed at",Sys.time(),"using", input$protein, "data"))
  })

}

shinyApp(ui = ui, server = server)