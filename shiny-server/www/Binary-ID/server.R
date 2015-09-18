##Loading the libraries
library(shiny)
library(poppr)
library(ape)
library(igraph)

#####################################################
# IMPORTANT: ALl the functions written before line 68 refer to functions loaded by R
# before shiny is deployed and executes the custom functions. These functions are definitions of
# global functions and variables. Make sure to add these kinds of functions here. For more information
# refer ro the shiny manual.
####################################################


########### MICROBE-ID customization  ############
# Here's where you add your database file (Comma Separated Object). Make sure
# that the database is in the same folder than this file (server.R)

df <- read.table("Aeut.txt", header = TRUE, sep = "\t")

# Transforming the dataframe into a matrix
df.m <- as.matrix(df)

########### MICROBE-ID customization  ############
# This line creates a function to obtain the genetic distances and add
# them in a switch (evaluates EXPR and accordingly chooses one of the further arguments )
# for easly post-hoc handling

get_dist_fun <- function(dist){
  switch(dist,
         nei = 'nei.dist',
         edwards = 'edwards.dist',
         rogers = 'rogers.dist',
         reynolds = 'reynolds.dist',
         provesti = 'provesti.dist')
}


# Functions to create elements to plot:

## 1. Distance Tree
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

## 2. Minimum spanning network
plot.minspan <- function(gen, mst, gadj=3, inds = "none", ...){
  plot_poppr_msn(gen, mst, gadj = gadj, vertex.label.color = "firebrick", inds = inds,
                 vertex.label.font = 2, vertex.label.dist = 0.5, nodelab = 100,
                 quantiles = FALSE)
}

########### MICROBE-ID customization  ############
# From this line on, every one of the functions is going to be used by shiny in a reactive way. All modifications
# of processes and outputs for the User interface file (in this case, the www/index.html) are found here.
#
# INPUT FROM index.html: All variables that start with index$...
# OUTPUT TO index.html: All variables that start with output$...
#
# To determine which variable communicates with whic <div> in the index.html file, search for the line with the
# class=shiny*.(e.g. The input$table variable is gonna be filled with info from the <div class="shiny-bound-input" id="table">.
#
# For more information refer to the shiny manual

shinyServer(function(input, output) {
  data.genoid <- reactive({
    if (gsub("\\s", "", input$table) == ""){
      return(NULL)
    } else {
      input_table <- read.table(text = input$table, stringsAsFactors = FALSE)
      if (input_table[1,1] == "Ind" | input_table[1,2] == "Pop"){
        input_table <- input_table[-1,]
      }
      colnames(input_table) <- colnames(df.m)
      input_data.genoid            <- input_table[[1]]
      df.m <- rbind(df.m, input_table, deparse.level = 0)
      df.m <- as.data.frame(df.m)
      gen  <- df2genind(df.m[, -c(1, 2)], ncode=1,  ploid = 2, pop = df.m[, 2], type="PA", ind.names = df.m[, 1])
    #Adding colors to the tip values according to the clonal lineage
      gen$other$tipcolor   <- pop(gen)
      gen$other$input_data.genoid <- input_data.genoid
      ngroups              <- length(levels(gen$other$tipcolor))
      ###########  MICROBE-ID customization ############
      # Change these colors to represent the groups defined in your data.genoid set.
      defined_groups <- c("blue", "darkolivegreen")
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

# Setting a random seed for the current session from the user interface (<input type = "number" name = "seed" id = "seed" value = "9449" min = "0" />)
  seed <- reactive({
    return(input$seed)
  })

# Matching the distances from the user interface (<select id="distance" style="width:300px">)
  distfun <- reactive({
    get_dist_fun(input$distance)
  })

# Greyscale slider settings from the user interface (<input id="integer" type="slider" name="integer" value="3" class="jslider" data-from="0" data-to="50" data-step="1" data-skin="plastic" data-round="FALSE" data-locale="us" data-format="#,##0.#####" data-smooth="FALSE"/>)
  slider <- reactive({
    slider.a <- (input$integer)
    return(slider.a)
  })

# Processing the results. The functions here create the figures to be displayed by the user interface.

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

############ MICROBE-ID customization ############
# The following lines of code communicate with the user interface to
# plot the outputs from the processes in the server script.

  ## Distance Tree 	(<div id="distPlotTree" class="span6 shiny-plot-output">)
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

  ##Minimum Spanning Network (<div id="MinSpanTree" class="shiny-plot-output")
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


############ MICROBE-ID customization ############
# The following lines of code communicate with the user interface to
# download the outputs from the processes in the server script.

  ## Distance tree in .tre format (<a id="downloadData" class="btn btn-primary shiny-download-link">")
  output$downloadData <- downloadHandler(
    filename = function() { paste0(input$tree, '.tre') },
    content = function(file) {
      write.tree(boottree(), file)
    })

  ## Distance tree in PDF format (	<a id="downloadPdf"  class="btn btn-info shiny-download-link">)
  output$downloadPdf <- downloadHandler(
    filename = function() { paste0(input$tree, '.pdf') },
    content = function(file) {
      pdf(file, width=11, height=8.5)
      plot.tree(tree, tip.col=as.character(unlist(gen$other)))
      dev.off()
    })

  ## Minimum spanning network in PDF format (<a id="downloadPdfMst"  class="btn btn-info shiny-download-link")
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
