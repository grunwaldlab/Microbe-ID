library(shiny)
library(poppr)
library(ape)
library(phangorn)
library(igraph)
library(gdata)

# Change this path to reflect where your binary of muscle is located.
#muscle_dir <- "/usr/bin/muscle"
muscle_dir <- "/Users/tabimaj/Downloads/muscle3.8.31_i86darwin64"

get_last_substring <- function(x, sep = "_"){
  splitx <- strsplit(x, sep)
  last   <- vapply(splitx, function(y) y[[length(y)]], character(1))
  return(last)
}

shinyServer(function(input, output) {
  
  


  # COMMENT:
  #
  # See www/index.html for my comments on this reactive block.
  # -Zhian

  
  ########### IMPORTANT ############
  # Here's where you add your database file (Comma Separated Object). Make sure 
  # that the database is in the same folder than this file (server.R)
  data_f <- reactive({ 
    df <- read.dna(input$dataset, format = "fasta")
    return(df)
    })

  alin <- reactive({
    if (gsub("\\s", "", input$fasta) == ""){
      return(NULL)
    } else {
      if (startsWith(input$fasta,">") == TRUE){
        cat(input$fasta, file="input.fasta")
        input_table <- read.dna("input.fasta", format="fasta")
        rownames(input_table) <- paste(rownames(input_table),c("query"),sep="_")
        input_table <- as.list(input_table)
        df <- as.list(data_f())
        all <- c(input_table,df)
        all.al <- muscle(all, exec = muscle_dir)
        return(all.al) # Returns a DNAbin object.
      }else{
        return(NULL)
      }
    }
    })
  

  dist <- reactive({
    all.dist <- dist.dna(alin(),input$model)
    return(all.dist)
  })
  
  data <- reactive({
    gen <- DNAbin2genind(alin())
    #Adding colors to the tip values according to the clonal lineage
    popnames <- get_last_substring(gen$ind.names, "_")
    pop(gen) <- popnames
    gen$other$tipcolor <- pop(gen)
    name_char_end <- nchar(indNames(gen)) - nchar(popnames) - 1
    gen$ind.names <- substr(indNames(gen), 1, name_char_end)
    gen$other$input_data <- indNames(gen)[pop(gen) %in% "query"]
    ngroups   <- length(levels(gen$other$tipcolor))
#     ########### IMPORTANT ############
#     # Change these colors to represent the groups defined in your data set.
#     #
     defined_groups <- c("blue", "darkcyan", "darkolivegreen", "darkgoldenrod")
#     #
#     # Change heat.colors to whatever color palette you want to represent
#     # submitted data. 
#     #
     input_colors   <- rainbow(ngroups - length(defined_groups))
#     #
#     ##################################
    levels(gen$other$tipcolor) <- c(defined_groups, input_colors)
     gen$other$tipcolor <- as.character(gen$other$tipcolor)
    return(gen)
  })
  
  seed <- reactive({
    return(input$seed)
  })
  
  boottree <- reactive({
    # Running the tree, setting a cutoff of 50 and saving it into a variable to 
    # be plotted (tree)
    if (input$boot > 1000){
      return(1000L)
    } else if (input$boot < 10){
      return(10L)
    }
    set.seed(seed())

    # COMMENT:
    # 
    # It looks like you had the right idea here on this commented line.
    # Remember, you can write functions that return functions.
    # - Zhian
  #  boot.dist <- function(al,di){
      if (input$tree == "upgma"){
        tre <- upgma(dist())
        bp <- boot.phylo(tre,alin(),function(x) upgma(dist.dna(x,input$model)),B=input$boot)
      } else {
        tre <- nj(dist())
        bp <- boot.phylo(tre,alin(),function(x) nj(dist.dna(x,input$model)),B=input$boot)
      }

      tre$node.labels <- round(((bp / input$boot)*100))
      #plot(tre,label.offset = 0.0125)
      #nodelabels(tre$node.labels, adj = c(1.3, -0.5), frame="n", cex=0.9, font=3)

    if (input$tree=="nj"){
      tre <- phangorn::midpoint(ladderize(tre))
    }
    return(tre)
  })
  
  msnet <- reactive ({
    msn.plot <-poppr.msn(data(), dist(), palette=rainbow, showplot = FALSE)
   # browser()
    V(msn.plot$graph)$size <- 3
    return(msn.plot)
  })


# COMMENT:
#
# Since we are drawing plots and saving those plots to files, they should be the
# same. Instead of copying/pasting the functions, we could write wrappers functions
# for these plotting functions Where the variables (data(), boottree(), msnet(), 
# seed(), etc.) are taken in as arguments to the function and options that we
# want to keep static (cex, border, etc.) are defined.
# - Zhian


  # COMMENT:
  #
  # Tiplabels are doubled up and not lining up with the tree. 
  # Remove the tiplabels command and use the tip.col argument in the
  # plot.phylo function. See the poppr internal function: poppr.plot.phylo
  # -Zhian

  output$distPlotTree <- renderPlot({
    if (is.null(alin())){
      plot.new() 
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) + 
      text(x = 0.5, y = 0.9, "No FASTA data has been input.", cex = 1.6, col = "white")
    } else if (is.integer(boottree())){
      msg <- ifelse(boottree() > 10L, "\nless than or equal to 1000",
                                      "greater than 10")
      msg <- paste("The number of bootstrap replicates should be", msg)
      plot.new()
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) + 
      text(x = 0.5, y = 0.9, msg, cex = 1.6, col = "white")
    } else {
        #Drawing the tree
      plot.phylo(boottree())
      
      #Adding the tip labels from each population, and with the already defined colors
      tiplabels(pop(data()), adj = c(-4, 0.5), frame = "n", 
                col = data()$other$tipcolor, cex = 0.8, font = 2)
      
      #Adding the nodel labels: Bootstrap values.
      nodelabels(boottree()$node.label, adj = c(1.2, -0.5), frame = "n", 
                 cex = 0.9, font = 3)
      
      if (input$tree == "upgma"){
        axisPhylo(3)
      } else {
        add.scale.bar(x = 0.89, y = 1.18, length = 0.05, lwd = 2)
      }
    }
  })
  
  #Minimum Spanning Network
  output$MinSpanTree <- renderPlot({
    if (is.null(alin())){
      plot.new() 
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) + 
      text(x = 0.5, y = 0.9, "No FASTA data has been input.", cex = 1.6, col = "white")
    } else {
      set.seed(seed())
      plot_poppr_msn(data(), msnet(), vertex.label.color = "firebrick", 
                     vertex.label.font = 2, vertex.label.dist = 0.5, 
 quantiles = FALSE, gadj = 40, inds = data()$other$input_data, nodelab = 10)
    }  	
    
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() { paste0(input$tree, '.tre') },
    content = function(file) {
      write.tree(boottree(), file)
    })
  
  output$downloadPdf <- downloadHandler(
    filename = function() { paste0(input$tree, '.pdf') },
    content = function(file) {
      pdf(file, width=11, height=8.5)
      plot.phylo(boottree(), cex = 0.5)
      tiplabels(pop(data()), adj = c(-4, 0.5), frame = "n", 
                col = data()$other$tipcolor, cex = 0.4, font = 2)
      nodelabels(boottree()$node.label, adj = c(1.2, -0.5), frame = "n", 
                 cex = 0.4, font = 3)
      if (input$tree == "upgma"){
        axisPhylo(3)
      }
      dev.off()
    })
  
  output$downloadPdfMst <- downloadHandler(
    filename = function() { paste0("min_span_net", '.pdf')} ,
    content = function(file) {
      pdf(file, width=11, height=8.5)
      set.seed(seed())
      plot_poppr_msn(data(), msnet(), vertex.label.color = "firebrick", 
                     vertex.label.font = 2, vertex.label.dist = 0.5, 
                     inds = data()$other$input_data, # This does not match with the above function.
                     quantiles = FALSE, gadj = 90)
      dev.off()
    }
  )
  
})
