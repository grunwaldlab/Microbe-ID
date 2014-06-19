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


# Functions to create elements to plot
## Distance Tree
plot.tree <- function (tree, type = input$tree, ...){
  ARGS <- c("nj", "upgma")
  type <- match.arg(type, ARGS)
  barlen <- min(median(tree$edge.length), 0.1)
  if (barlen < 0.1) 
    barlen <- 0.01
  plot.phylo(tree, cex = 0.8, font = 2, adj = 0, xpd = TRUE, 
             label.offset = 0.0125, ...)
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", 
             cex = 0.8, font = 3, xpd = TRUE)
  if (type == "nj") {
    add.scale.bar(lwd = 5, length = barlen)
    tree <- ladderize(tree)
  }
  else {
    axisPhylo(3)
  }
}
## Minimum spanning network
plot.minspan <- function(gen, mst, gadj=3, inds = "none", ...){
  plot_poppr_msn(gen, mst, gadj = gadj, vertex.label.color = "firebrick", inds = inds,
                 vertex.label.font = 2, vertex.label.dist = 0.5, nodelab = 100,
                 quantiles = FALSE)  
}

shinyServer(function(input, output) {
  

  data.genoid <- reactive({
    if (gsub("\\s", "", input$table) == ""){
      return(NULL)
    } else {
      #browser()
      input_table <- read.table(text = input$table, stringsAsFactors = FALSE)
      if (input_table[1,1] == "Ind" | input_table[1,2] == "Pop"){
        input_table <- input_table[-1,]
      }
      colnames(input_table) <- colnames(df.m)
      input_data.genoid            <- input_table[[1]]
      df.m <- rbind(df.m, input_table, deparse.level = 0)
      df.m <- as.data.frame(df.m)
      gen  <- df2genind(df.m[, -c(1, 2)], ploid = 2, pop = df.m[, 2], type='PA', ind.names = df.m[, 1], )
      #Adding colors to the tip values according to the clonal lineage
      gen$other$tipcolor   <- pop(gen)
      gen$other$input_data.genoid <- input_data.genoid
      ngroups              <- length(levels(gen$other$tipcolor))
      ########### IMPORTANT ############
      # Change these colors to represent the groups defined in your data.genoid set.
      #
      defined_groups <- c("blue", "darkolivegreen")#, "red")
      #
      # Change heat.colors to whatever color palette you want to represent
      # submitted data.genoid. 
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
 
# Random Seed
  seed <- reactive({
    return(input$seed)
  })
   
# Mathcing the distances from the UI to the server end
  distfun <- reactive({
    get_dist_fun(input$distance)
  })

# Greyscale slider
  slider <- reactive({
    slider.a <- (input$integer)
    return(slider.a)
  })

# Calculating the results

  # Distance tree with bootstrap
  boottree <- reactive({
      if (input$boot > 1000){
        return(1000L)
      } else if (input$boot < 10){
        return(10L)
      }
      set.seed(seed())      
      DIST <- match.fun(distfun())
      tree <- aboot(data.genoid(), distance = DIST, sample = input$boot, showtree=FALSE,  
                    tree = input$tree, cutoff = 50)
      if (input$tree=="nj"){
        tree <- phangorn::midpoint(ladderize(tree))
      }
      tree$tip.label <- paste(tree$tip.label, as.character(pop(data.genoid()))) 
      return(tree)
    })
  
  # Minimum spanning network
    msnet <- reactive ({
      DIST <- match.fun(distfun())
      msn.plot <- poppr.msn(data.genoid(),distmat=DIST(data.genoid()),showplot=FALSE)
      V(msn.plot$graph)$size <- 3
      return(msn.plot)
    })

# Plotting on the UI
  
  ## Distance Tree
  output$distPlotTree <- renderPlot({
    if (is.null(data.genoid())){
      plot.new() 
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) + 
        text(x = 0.5, y = 0.9, "No data genoid has been input.", cex = 1.6, col = "white")
    } else if (is.integer(boottree())){
      msg <- ifelse(boottree() > 10L, "\nless than or equal to 1000",
                    "greater than 10")
      msg <- paste("The number of bootstrap replicates should be", msg)
      plot.new()
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) + 
        text(x = 0.5, y = 0.9, msg, cex = 1.6, col = "white")
    } else {
      plot.tree(boottree(), type = input$tree, tip.col=as.character(unlist(data.genoid()$other$tipcolor)))
    }
  })
  
  ##Minimum Spanning Network
  output$MinSpanTree <- renderPlot({
    if (is.null(data.genoid())){
      plot.new() 
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) + 
        text(x = 0.5, y = 0.9, "No data has been input.", cex = 1.6, col = "white")
    } else {
      set.seed(seed())
      plot.minspan(data.genoid(), msnet(), gadj = slider(), inds = data.genoid()$other$input_data.genoid)
    }
  })
  
  
#Downloading results
  
  ## Distance tree in .tre format
  output$downloadData <- downloadHandler(
    filename = function() { paste0(input$tree, '.tre') },
    content = function(file) {
      write.tree(boottree(), file)
    })
  
  ## Distance tree in PDF format  
  output$downloadPdf <- downloadHandler(
    filename = function() { paste0(input$tree, '.pdf') },
    content = function(file) {
      pdf(file, width=11, height=8.5)
      plot.tree(tree, tip.col=as.character(unlist(gen$other)))
      dev.off()
    })
  
  ## Minimum spanning network in PDF format
  output$downloadPdfMst <- downloadHandler(
    filename = function() { paste0("min_span_net", '.pdf')} ,
    content = function(file) {
      pdf(file, width=11, height=8.5)
      set.seed(seed())
      plot.minspan(data.genoid(),msnet())
      dev.off()
    }
  )
  
  #EOF
})

