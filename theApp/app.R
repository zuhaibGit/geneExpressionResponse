library(shiny)
library(ggplot2)
library(stringr)
library(ggplot2)
library(BiocManager)
#options(repos = BiocManager::repositories())
setwd("/home/zuhaib/Desktop/covid19Research/geneResponseApp/geneExpressionResponse/theApp")
# Reads in the long_ files
fls <- unlist(lapply(list.files()[grep("data_", list.files())], function(d) {
  path <- paste0("./", d, "/")
  filesInDir <- list.files(path)
  return(paste0(path, filesInDir[grep("long_", filesInDir)]))
}))
datasets <- lapply(fls, function(x) {
  read.table(x, header = T, sep = "\t")
})
names(datasets) <- fls
names(datasets) <- str_replace_all(names(datasets), "\\.txt", "")
names(datasets) <- str_replace_all(names(datasets), "\\..+long_", "")
# Need to remove the decimal part of the GSE148729 dataset gene IDs, to make it easier to search
datasets$GSE148729_Calu3_mockInfection$Gene <-  str_replace_all(datasets$GSE148729_Calu3_mockInfection$Gene, "\\..*", "")
datasets$GSE148729_Calu3_sarsCov1$Gene <- str_replace_all(datasets$GSE148729_Calu3_sarsCov1$Gene, "\\..*", "")
datasets$GSE148729_Calu3_sarsCov2$Gene <- str_replace_all(datasets$GSE148729_Calu3_sarsCov2$Gene, "\\..*", "")


# Reads in the normalized expression data
fls2 <- unlist(lapply(list.files()[grep("data_", list.files())], function(d) {
  path <- paste0("./", d, "/")
  filesInDir <- list.files(path)
  return(paste0(path, filesInDir[grep("normalized_", filesInDir)]))
}))
normData <- lapply(fls2, function(x) {
  read.table(x, header = T, sep = "\t")
})
names(normData) <- fls2
names(normData) <- str_replace_all(names(normData), "\\.txt", "")
names(normData) <- str_replace_all(names(normData), "\\..+normalized_", "")

fls3 <- paste0("./geneMappings/", list.files("./geneMappings")[grep("mappings_", list.files("./geneMappings"))])
mappings <- lapply(fls3, function(x) {
  read.table(x, header = T, sep = "\t")
})
names(mappings) <- fls3
names(mappings) <- str_replace_all(names(mappings), "\\.txt", "")
names(mappings) <- str_replace_all(names(mappings), "\\..+mappings_", "")

# Reads in file specifying which samples correspond to which columns in the normalized data.
# This is the help disply the desired data in tabular format
groupToColumns <- read.table("groupToColumns.txt", header = T, sep = "\t")
groupToColumns$File <- str_replace_all(groupToColumns$File, "\\.txt", "")
groupToColumns$File <- str_replace_all(groupToColumns$File, ".+normalized_", "")


# adf <- datasets$GSE148729_Calu3_sarsCov2
# bdf <- mappings$ensembl
# testEntrez <- 6013
# g <- bdf[which(bdf$ENTREZ_ID == testEntrez),1]
# adf[grep(g[1], adf[,1]),]
# adf[grep(g[2], adf[,1]),]

# Define UI for miles per gallon app ----
ui <- pageWithSidebar(
  
  # App title ----
  headerPanel("Gene Expression Changes"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    checkboxGroupInput("selDatasets", "Select Datasets", choices = names(datasets)),
    radioButtons("selID", label = "What type of gene ID?", 
                 choices = list("Gene Symbol" = "GENE_SYMBOL", "Ensembl ID" = "ENSEMBL_GENE_ID", "Entrez ID" = "ENTREZ_ID")),
    textAreaInput("selGenes", "Genes of Interest", height = "200px"),
    radioButtons("collapseLines", label = "Collapse lines?",
                 choices = list("Yes" = "yes", "No" = "no"), 
                 selected = "no"),
    submitButton("Submit")
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    verbatimTextOutput("GNF"),
    tabsetPanel(type = "tabs", 
                tabPanel("main", uiOutput("Plots")),
                tabPanel("table", uiOutput("Tables")))
  )
)

# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  ########## FUNCTIONS ##########
  # Takes in the DE between time points of some gene, and returns x,y coordinates for the line
  # as well as the color of the points based on whether it was significantly expressed.
  # Note: Time points must be sorted
  makeLine <- function(timePoints, collapseLinesFlag = F) {
    minTimePoint <- timePoints[1,2]
    # If the lines will collapse, then raise each y value to the power of to so lines don't clump together
    if (collapseLinesFlag == T) {
      yVals <- c(1, cumprod(2^(timePoints$log2FoldChange)))
    } else {
      yVals <- c(0, cumsum(timePoints$log2FoldChange))
    }
    retDF <- data.frame(x = c(minTimePoint, timePoints$T2),
                        y = yVals,
                        sig = c("Significant", timePoints$Colour),
                        gene = timePoints[1,1])
    return(retDF)
  }
  ########## FUNCTIONS ##########
  
  
  # Based on the user-selected genes and datasets, creates a line for each gene in each dataset
  # Returns the line for each gene in each dataset (lns), also returns the max and min x and y values (useful for defining plot dimentions)
  # as well as the name of the datset (Name)
  dataToPlot <- reactive({
    if (input$collapseLines == "yes") {
      collapseFlag <- 0
    } else {
      collapseFlag <- 1
    }
    return(lapply(input$selDatasets, function(y) {
      ds <- datasets[[y]]
      genesOfInterest <- strsplit(input$selGenes, split = '[\r\n]')[[1]] 
      selectedIDType <- input$selID
      dataset_primaryID <- groupToColumns$PrimaryID[which(groupToColumns$Group == y)]
      lst <- lapply(genesOfInterest, function(x) {
        # First, map the users query genes to the gene IDs present in their dataset of interest
        m <- mappings[[dataset_primaryID]]
        g <- unique(m[grep(paste0("^", x, "$"), m[[selectedIDType]]),1])
        retLst <- lapply(g, function(z) {
          ret <- ds[grep(paste0("^", z, "$"), ds[,1]),]
          if (nrow(ret) == 0) {
            return(NA)
          } else {
            ret[,1] <- x
            return(ret)
          }
        })
        retLst <- retLst[which(!(is.na(retLst)))]
        #print(retLst)
        if (length(retLst) == 0) {
          return(NA)
        } else {
          return(retLst)
        }
      })
      genesNotFound <- genesOfInterest[which(is.na(lst))]
      lst <- lst[which(!is.na(lst))]
      lst <- unlist(lst, recursive = F)
      lns <- lapply(lst, function(x) {
        if (input$collapseLines == "yes") {
          return(makeLine(x, T))
        } else {
          return(makeLine(x, F))
        }
      })
      xMax <- max(unlist(lapply(lns, function(x) return(x[,1])))) + 1
      yMin <- min(unlist(lapply(lns, function(x) return(x[,2])))) - 1
      yMax <- max(unlist(lapply(lns, function(x) return(x[,2])))) + length(lns)*collapseFlag
      return(list(Plot = lns, xMax = xMax, yMin = yMin, yMax = yMax, Name = y, notFound = genesNotFound, DS = y))
    }))
  })
  
  # Defines the plotting windows based on the max and min x and y values for all plots
  # This is so the plot scales are consistent.
  plotWindow <- reactive({
    maxX <- max(unlist(lapply(dataToPlot(), function(x) return(x$xMax))))
    minX <- -9
    minY <- min(unlist(lapply(dataToPlot(), function(x) return(x$yMin))))
    maxY <- max(unlist(lapply(dataToPlot(), function(x) return(x$yMax))))
    if (input$collapseLines == "yes") {
      maxX <- maxX + 6
      minX <- -0.5
    } 
    
    return(list(right = maxX, left = minX, top = maxY, bottom = minY))
  })
  
  output$GNF <- renderPrint({
    lapply(dataToPlot(), function(d) {
      return(paste("In", d$DS, "we didn't find", d$notFound))
    })
  })
  
  # Renders the plots. One plot for each selected dataset. The plot scales based on the max and min y-values
  output$Plots <- renderUI({
      lapply(dataToPlot(), function(d) {
        if (input$collapseLines == "yes") {
          collapseFlag <- 0
        } else {
          collapseFlag <- 1
        }
        
        vShift <- 0
        renderPlot({
          plot(1, type="n", xlab="Time (h)", ylab="", xlim=c(plotWindow()$left, plotWindow()$right), ylim=c(plotWindow()$bottom, plotWindow()$top), main = d$Name)
          for (i in d$Plot) {
            lines(i$x,
                  i$y + vShift,
                  type = "b",
                  col = sapply(i$sig, function(x) if (x == "Significant") return("black") else return("yellow")),
                  cex = 2,
                  pch = 16)
            text(-5*collapseFlag + (1 - collapseFlag)*((i$x[length(i$x)]) + 3), 
                 vShift*collapseFlag + (1 - collapseFlag)*(i$y[length(i$y)]), 
                 labels = i$gene)
            vShift <- vShift + collapseFlag
          }
        }, height = max(20, (plotWindow()$top - plotWindow()$bottom) * (3100/(plotWindow()$top - plotWindow()$bottom) + 800) / 35))
      })
  })
  
  output$Tables <- renderUI({
    lapply(dataToPlot(), function(d) {
      renderTable({
        selDataSet <- d$Name
        selNormData <- groupToColumns[which(groupToColumns$Group == selDataSet),2]
        selCols <- c(1,as.integer(strsplit(groupToColumns[which(groupToColumns$Group == selDataSet),3], ",")[[1]]))
        retDF <- normData[[selNormData]]
        return(retDF[1:3,selCols])
      })
    })
  })
}

shinyApp(ui, server)









