library(shiny)
library(poppr)
library(ape)
library(igraph)

########### IMPORTANT ############
# Here's where you add your database file (Comma Separated Object). Make sure 
# that the database is in the same folder than this file (server.R)
df <- read.table("Aeut.txt", header = TRUE, sep = "\t")
##################################
df.m <- as.matrix(df)


get_dist_fun <- function(dist){
  switch(dist,
         nei = 'nei.dist',
         edwards = 'edwards.dist',
         rogers = 'rogers.dist',
         reynolds = 'reynolds.dist',
         provesti = 'provesti.dist')
}

shinyServer(function(input, output) {
  

  data <- reactive({
    if (gsub("\\s", "", input$table) == ""){
      return(NULL)
    } else {
      #browser()
      input_table <- read.table(text = input$table, stringsAsFactors = FALSE)
      if (input_table[1,1] == "Ind" | input_table[1,2] == "Pop"){
        input_table <- input_table[-1,]
      }
      colnames(input_table) <- colnames(df.m)
      input_data            <- input_table[[1]]
      df.m <- rbind(df.m, input_table, deparse.level = 0)
      df.m <- as.data.frame(df.m)
      gen  <- df2genind(df.m[, -c(1, 2)], ploid = 2, pop = df.m[, 2], type='PA', ind.names = df.m[, 1], )
      #Adding colors to the tip values according to the clonal lineage
      gen$other$tipcolor   <- pop(gen)
      gen$other$input_data <- input_data
      ngroups              <- length(levels(gen$other$tipcolor))
      ########### IMPORTANT ############
      # Change these colors to represent the groups defined in your data set.
      #
      defined_groups <- c("blue", "darkolivegreen")#, "red")
      #
      # Change heat.colors to whatever color palette you want to represent
      # submitted data. 
      #
      input_colors   <- heat.colors(ngroups - length(defined_groups))
      #
      ##################################
      levels(gen$other$tipcolor) <- c(defined_groups, input_colors)
      gen$other$tipcolor <- as.character(gen$other$tipcolor)
      gen <- as.genclone(gen)
      return(gen)
    }
  })
  
  seed <- reactive({
    return(input$seed)
  })
  
 

  distfun <- reactive({
    get_dist_fun(input$distance)
  })
  

boottree <- reactive({
    if (input$boot > 1000){
      return(1000L)
    } else if (input$boot < 10){
      return(10L)
    }
    set.seed(seed())
    
    DIST <- match.fun(distfun())
    tree <- aboot(data(), distance = DIST, sample = input$boot, showtree=FALSE,  
                  tree = input$tree, cutoff = 50)
    # This is a catch to avoid having missing data within the distance matrix. 
#     if ("try-error" %in% class(tree)){
#       for (i in sample(100)){
#         browser()
#         tree <- try(aboot(data(), dist_mat(), sample = input$boot, 
#                           tree = input$tree, cutoff = 50), silent = TRUE)
#         if (!"try-error" %in% class(tree)){
#           print(paste0("Success: ", i))
#           break
#         }
#         print(paste0("Failed: ", i))
#       }
    #}

    if (input$tree=="nj"){
      tree <- phangorn::midpoint(ladderize(tree))
    }
    return(tree)
  })
  
slider <- reactive({
  slider.a <- (input$integer)
  return(slider.a)
})

  msnet <- reactive ({
    DIST <- match.fun(distfun())
    msn.plot <- poppr.msn(data(),distmat=DIST(data()),showplot=FALSE)
    V(msn.plot$graph)$size <- 3
    return(msn.plot)
  })

  
  output$distPlotTree <- renderPlot({
    if (is.null(data())){
      plot.new() 
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) + 
      text(x = 0.5, y = 0.9, "No Binary (AFLP) data has been input.", cex = 1.6, col = "white")
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
      tiplabels(pop(data()), adj = c(-1, 0.5), frame = "n", 
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
    if (is.null(data())){
      plot.new() 
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) + 
      text(x = 0.5, y = 0.9, "No binary (AFLP) data has been input.", cex = 1.6, col = "white")
    } else {
      set.seed(seed())
      plot_poppr_msn(data(), msnet(), gadj=c(slider()), vertex.label.color = "firebrick", 
                     vertex.label.font = 2, vertex.label.dist = 0.5, 
                     inds = data()$other$input_data, quantiles = FALSE, nodelab = 10)
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
                     inds = data()$other$input_data, quantiles = FALSE)
      dev.off()
    }
  )
  
})
