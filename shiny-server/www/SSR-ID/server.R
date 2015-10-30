library(shiny)
library(poppr)
library(phangorn)
library(ape)
library(igraph)

#####################################################
# IMPORTANT: ALl the functions written before line 68 refer to functions loaded by R
# before shiny is deployed and executes the custom functions. These functions are definitions of'
# global functions and variables. Make sure to add these kinds of functions here. For more information
# refer ro the shiny manual.
####################################################


########### MICROBE-ID customization ############
# Here's where you add your database file (Comma Separated Object). Make sure
# that the database is in the same folder than this file (server.R)
df <- read.table("Ramorum_ssr.csv", header = TRUE, sep = "\t")
##################################
df.m <- as.matrix(df)



########### MICROBE-ID customization ############
# Change these values to the repeat lenghts and names of your SSR markers.
ssr <- c(PrMS6       = 3,
         PRMS9c3     = 2,
         PrMS39a     = 2,
         PrMS45      = 4,
         KI18        = 2,
         KI64        = 2,
         ILVOPrMS131 = 2
         )
##################################

# Functions to create elements to plot
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
  plot_poppr_msn(gen, mst, gadj=gadj, vertex.label.color = "firebrick", inds = inds,
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

#Data input and manipulation
  data.genoid <- reactive({
    if (gsub("\\s", "", input$table) == ""){
      return(NULL)
    } else {
      input_table <- read.table(text = input$table, stringsAsFactors = FALSE)
      colnames(input_table) <- colnames(df.m)
      input_data.genoid            <- input_table[[1]]
      df.m <- rbind(df.m, input_table, deparse.level = 0)
      df.m <- as.data.frame(df.m)
      gen  <- df2genind(df.m[, -c(1, 2)], ploid = 2, sep = "/", pop = df.m[, 2],
                        ind.names = df.m[, 1])
      #Adding colors to the tip values according to the clonal lineage
      gen$other$tipcolor   <- pop(gen)
      gen$other$input_data.genoid <- input_data.genoid
      ngroups              <- length(levels(gen$other$tipcolor))
      ########### MICROBE-ID customization ############
      # Change these colors to represent the groups defined in your data.genoid set.
      #
      defined_groups <- c("blue", "darkcyan", "darkolivegreen", "darkgoldenrod","red")
      #
      # Change heat.colors to whatever color palette you want to represent
      # submitted data.genoid.
      #
      input_colors   <- heat.colors(ngroups - length(defined_groups))
      #
      ##################################
      levels(gen$other$tipcolor) <- c(defined_groups, input_colors)
      gen$other$tipcolor <- as.character(gen$other$tipcolor)
      return(gen)
    }
  })

# Setting a random seed for the current session from the user interface (<input type = "number" name = "seed" id = "seed" value = "9449" min = "0" />)
  seed <- reactive({
    return(input$seed)
  })

# Greyscale slider settings from the user interface (<input id="integer" type="slider" name="integer" value="3" class="jslider" data-from="0" data-to="50" data-step="1" data-skin="plastic" data-round="FALSE" data-locale="us" data-format="#,##0.#####" data-smooth="FALSE"/>)

  slider <- reactive({
    slider.a <- (input$integer)
    return(slider.a)
  })

# Processing the results. The functions here create the figures to be displayed by the user interface.
# Bootstrap of a distance tree out of the data
  boottree <- reactive({
    # Running the tree, setting a cutoff of 50 and saving it into a variable to
    # be plotted (tree)
    if (input$boot > 1000){
      return(1000L)
    } else if (input$boot < 10){
      return(10L)
    }
    set.seed(seed())
    tree <- try(bruvo.boot(data.genoid(), replen = ssr, sample = input$boot,
                       tree = input$tree, cutoff = 50), silent = TRUE)

    # This is a catch to avoid having missing data within the distance matrix.
    if ("try-error" %in% class(tree)){
      for (i in sample(100)){
        tree <- try(bruvo.boot(data.genoid(), replen = ssr, sample = input$boot,
                               tree = input$tree, cutoff = 50), silent = TRUE)
        if (!"try-error" %in% class(tree)){
          print(paste0("Success: ", i))
          break
        }
        print(paste0("Failed: ", i))
      }
    }
    tree$tip.label <- paste(tree$tip.label, as.character(pop(data.genoid())))
    if (input$tree=="nj"){
      tree <- phangorn::midpoint(ladderize(tree))
    }
    return(tree)
  })

#Minimum spanning network creation
  msnet <- reactive ({
    msn.plot <- bruvo.msn(data.genoid(), replen = ssr)
    V(msn.plot$graph)$size <- 3
    return(msn.plot)
  })


############ MICROBE-ID customization ############
# The following lines of code communicate with the user interface to
# plot the outputs from the processes in the server script.

## Distance Tree (<div id="distPlotTree" class="span6 shiny-plot-output">)
  output$distPlotTree <- renderPlot({
    if (is.null(data.genoid())){
      plot.new()
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) +
      text(x = 0.5, y = 0.9, "No SSR data has been input.", cex = 1.6, col = "white")
    } else if (is.integer(boottree())){
      msg <- ifelse(boottree() > 10L, "\nless than or equal to 1000",
                                      "greater than 10")
      msg <- paste("The number of bootstrap replicates should be", msg)
      plot.new()
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) +
      text(x = 0.5, y = 0.9, msg, cex = 1.6, col = "white")
    } else {
      plot.tree(boottree(), type=input$tree, tip.col=as.character(unlist(data.genoid()$other$tipcolor)))
    }
  })

##Minimum Spanning Network (<div id="MinSpanTree" class="shiny-plot-output")
  output$MinSpanTree <- renderPlot({
    if (is.null(data.genoid())){
      plot.new()
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) +
      text(x = 0.5, y = 0.9, "No SSR data has been input.", cex = 1.6, col = "white")
    } else {
      set.seed(seed())
      plot.minspan(data.genoid(), msnet(), gadj=slider(), inds = data.genoid()$other$input_data.genoid)
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

## Distance tree in PDF format (<a id="downloadPdf"  class="btn btn-info shiny-download-link">)
  output$downloadPdf <- downloadHandler(
    filename = function() { paste0(input$tree, '.pdf') },
    content = function(file) {
      pdf(file, width=11, height=8.5)
      plot.tree(boottree(), type=input$tree, tip.col=as.character(unlist(data.genoid()$other$tipcolor)))
      dev.off()
    })

## Minimum spanning network in PDF format (<a id="downloadPdfMst"  class="btn btn-info shiny-download-link">)
  output$downloadPdfMst <- downloadHandler(
    filename = function() { paste0("min_span_net", '.pdf')} ,
    content = function(file) {
      pdf(file, width=11, height=8.5)
      set.seed(seed())
      plot.minspan( data.genoid(), msnet(), gadj=slider(), inds = data.genoid()$other$input_data.genoid)
      dev.off()
    }
  )

#EOF
})
