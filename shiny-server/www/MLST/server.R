library(shiny)
library(poppr)
library(ape)
library(phangorn)
library(igraph)
library(gdata)

########### IMPORTANT ############
# Here's where you add your database file (Comma Separated Object). Make sure 
# that the database is in the same folder than this file (server.R)
df <- read.dna("database.fasta", format="fasta")




#==============================================================================#
# Function to plot the minimum spanning network obtained from poppr.msn or
# bruvo.msn. This is in development. At the moment, it's purpose is to plot all
# of the individual names that are assigned to a node. 
# usage:
# plot_poppr_msn(gid, poppr_msn)
# gid: a genind object
# poppr_msn: the output from a poppr.msn or bruvo.msn function
# gadj: the grey adjustment factor
# glim = grey limit
# gweight = -1 for greater distances, 1 for shorter distances
# inds = a character vector containing names of individuals you want to label on the graph.
#        set to "none" or any character that doesn't exist in your data for no labels
# quantiles = TRUE if you want the grey scale to be plotted using the raw data. 
#             FALSE if you want the range from minimum to maximum plotted.
# nodelab = when there are no inds labeled, this will label nodes with the number of individuals
#           greater than or equal to this number. 
# 
# Example:
#
# library(poppr)
# source('poppr_plot_msn.r')
# data(partial_clone)
# pc.msn <- bruvo.msn(partial_clone, replen = rep(1, 10))
# plot_poppr_msn(partial_clone, pc.msn, vertex.label.color = "firebrick", vertex.label.font = 2)
#==============================================================================#

plot_poppr_msn <- function(gid, poppr_msn, gadj = 3, glim = c(0, 0.8),
                           gweight = 1, inds = "ALL", quantiles = TRUE, nodelab = 2, ...){
  if(!is.genind(gid)) stop(paste(gid, "is not a genind object."))
  if(!identical(names(poppr_msn), c("graph", "populations", "colors"))) stop("graph not compatible")
  
  # Importing functions from igraph. This can be deleted when implemented in
  # poppr.
  E <- match.fun(igraph::E)
  `E<-` <- match.fun(igraph::`E<-`)
  V <- match.fun(igraph::V)
  plot.igraph <- match.fun(igraph::plot.igraph)
  
  # Adjusting color scales. This will replace any previous scaling contained in
  # poppr_msn.
  weights <- E(poppr_msn$graph)$weight
  wmin <- min(weights)
  wmax <- max(weights)
  gadj <- ifelse(gweight == 1, gadj, -gadj)
  E(poppr_msn$graph)$color <- gray(poppr:::adjustcurve(weights, glim=glim,
                                                       correction=gadj,show=FALSE))
  
  # Highlighting only the names of the submitted genotypes and the matching
  # isolates.
  gid.mlg <- mlg.vector(gid)
  labs <- unique(gid.mlg)
  # The labels in the graph are organized by MLG, so we will use that to extract
  # the names we need.
  if (length(inds) == 1 & toupper(inds[1]) == "ALL"){
    gid.input <- unique(gid.mlg)
  } else {
    gid.input <- unique(gid.mlg[gid@ind.names %in% inds])
  }
  # Combine all the names that match with each particular MLG in gid.input.
  combined_names <- vapply(gid.input, function(x)
    paste(rev(gid@ind.names[gid.mlg == x]),
          collapse = "\n"),
    " ")
  # Remove labels that are not specified.
  labs[which(!labs %in% gid.input)] <- NA
  labs[!is.na(labs)] <- combined_names
  if (all(is.na(labs))){
    labs <- V(poppr_msn$graph)$size
    labs <- ifelse(labs >= nodelab, labs, NA)
  }
  # Change the size of the vertices to a log scale.
  vsize <- log(V(poppr_msn$graph)$size, base = 1.15) + 3
  
  # Plotting parameters.
  def.par <- par(no.readonly = TRUE)
  # Setting up the matrix for plotting. Three vertical panels with a 1:3:1
  # ratio with legend:plot:greyscale
  #layout(matrix(1:3, ncol = 3), widths = c(1, 3, 0.5))
  layout(matrix(c(1,2,1,3), ncol = 2, byrow = TRUE),
         widths = c(1, 4), heights= c(4.5, 0.5))
  # mar = bottom left top right
  
  ## LEGEND
  par(mar = c(0, 0, 1, 0) + 0.5)
  too_many_pops <- as.integer(ceiling(length(gid$pop.names)/30))
  pops_correction <- ifelse(too_many_pops > 1, -1, 1)
  yintersperse <- ifelse(too_many_pops > 1, 0.51, 0.62)
  plot(c(0, 2), c(0, 1), type = 'n', axes = F, xlab = '', ylab = '',
       main = 'POPULATION')
  legend("topleft", bty = "n", cex = 1.2^pops_correction,
         legend = poppr_msn$populations, fill = poppr_msn$color, border = NULL,
         ncol = too_many_pops, x.intersp = 0.45, y.intersp = yintersperse)
  
  ## PLOT
  par(mar = c(0,0,0,0))
  plot.igraph(poppr_msn$graph, vertex.label = labs, vertex.size = vsize, ...)
  
  ## SCALE BAR
  if (quantiles){
    scales <- sort(weights)
    greyscales <- gray(poppr:::adjustcurve(scales, show=FALSE,
                                           glim=glim, correction=gadj))
  } else {
    scales <- seq(wmin, wmax, l = 1000)
    greyscales <- gray(poppr:::adjustcurve(scales, show=FALSE,
                                           glim=glim, correction=gadj))
  }
  legend_image <- as.raster(matrix(greyscales, nrow=1))
  par(mar = c(0, 1, 0, 1) + 0.5)
  plot.new()
  rasterImage(legend_image, 0, 0.5, 1, 1)
  polygon(c(0, 1 , 1), c(0.5, 0.5, 0.8), col = "white", border = "white", lwd = 2)
  axis(3, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(scales), 3))
  text(0.5, 0, labels = "DISTANCE", font = 2, cex = 1.5, adj = c(0.5, 0))
  
  # Return top level plot to defaults.
  layout(matrix(c(1), ncol=1, byrow=T))
  par(mar=c(5,4,4,2) + 0.1) # number of lines of margin specified.
  par(oma=c(0,0,0,0)) # Figure margins
}


shinyServer(function(input, output) {
  

  alin <- reactive({
    if (gsub("\\s", "", input$fasta) == ""){
      return(NULL)
    } else {
      if (startsWith(input$fasta,">") == TRUE){
        cat(input$fasta,file="input.fasta")
        input_table <- read.dna("input.fasta",format="fasta")
        rownames(input_table) <- paste(rownames(input_table),c("query"),sep="_")
        input_table <- as.list(input_table)
        df <-as.list(df)
        all <- c(input_table,df)
        all.al <- muscle(all,exec="/Users/tabimaj/Downloads/muscle3.8.31_i86darwin64")
        return(all.al)
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
    pop(gen) <- as.factor(sapply(strsplit(gen$ind.names,"_"),"[[",2))
    gen$other$tipcolor <- pop(gen)
    gen$ind.names <- sapply(strsplit(gen$ind.names,"_"),"[[",1)
    gen$other$input_data <- alin()
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
  #  boot.dist <- function(al,di){
      if(input$tree == "upgma"){
      tre <- upgma(dist())
      bp <- boot.phylo(tre,alin(),function(x) upgma(dist()),B=input$boot)
      }else{
      tre <- nj(dist())
      bp <- boot.phylo(tre,alin(),function(x) nj(dist()),B=input$boot)
      }

      tre$node.labels <- round(((bp / input$boot)*100))
      plot(tre,label.offset = 0.0125)
      nodelabels(tre$node.labels, adj = c(1.3, -0.5), frame="n", cex=0.9, font=3)

    if (input$tree=="nj"){
      tre <- phangorn::midpoint(ladderize(tre))
    }
    return(tre)
  })
  
  msnet <- reactive ({
    msn.plot <-poppr.msn(data(), dist(),palette=rainbow)
   # browser()
    V(msn.plot$graph)$size <- 3
    return(msn.plot)
  })

  
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
 quantiles = FALSE, gadj = 40)
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
                     inds = data()$other$input_data, quantiles = FALSE, gadj = 90)
      dev.off()
    }
  )
  
})
