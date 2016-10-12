library(shiny)
library(poppr)
library(ape)
library(phangorn)
library(igraph)
library(gdata)
library(XML , lib.loc="/data/www/chang/r_libs/")
library(phyloch, lib.loc="/data/www/chang/r_libs/")

#####################################################
# IMPORTANT: ALl the functions written before line 69 refer to functions loaded by R
# before shiny is deployed and executes the custom functions. These functions are definitions of
# global functions and variables. Make sure to add these kinds of functions here. For more information
# refer ro the shiny manual.
####################################################

########### MICROBE-ID customization  ############
# Change this path to reflect where your binary of mafft is located.
mafft_bin <- "/usr/bin/mafft"

########### MICROBE-ID customization  ############
#The name of the user submitted isolate for labelling in the tree
query_name <- "query_isolate"

########### MICROBE-ID customization  ############
#The dataset gene names. Add more for each new dataset. Each gene name refers to a file matching that name, ie GyrB.aln
#Alignments for each dataset are in folders named by that dataset, ie the index.html dataset option test-dataset refers to
#a folder containing the gene alignments for that dataset.

test_dataset_defined_genes <- c("CelA", "DnaA", "GyrB", "KdpA", "LigA", "NagA", "SdhA", "TomA")


get_last_substring <- function(x, sep = "_"){
  splitx <- strsplit(x, sep)
  last   <- vapply(splitx, function(y) y[[length(y)]], character(1))
  return(last)
}



# Functions to create elements to plot
## 1. Plotting trees
plot.tree <- function (tree, type = input$tree, ...){
  ARGS <- c("nj", "njs", "upgma")
  type <- match.arg(type, ARGS)
  barlen <- min(median(tree$edge.length), 0.1)
  if (barlen < 0.1)
    barlen <- 0.01
  plot.phylo(tree, type = "phylogram", cex = 0.8, font = 2, adj = 0, xpd = TRUE, label.offset = 0, ...)
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

## Plotting Minimum spanning network
plot.minspan <- function(gen, mst, gadj=3, inds="none", ...){
  plot_poppr_msn(gen, mst, gadj=gadj, vertex.label.color = "firebrick",nodelab=100,
                 vertex.label.font = 2, vertex.label.dist = 0.5,
                 inds = inds, quantiles = FALSE)
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

shinyServer(function(input, output, session) {

  data_f <- reactive({
    df <- read.dna(input$dataset, format = "fasta")
    return(df)
    })

alin <- reactive({
    if (gsub("\\s", "", input$fasta) == ""){
      return(NULL)
    } else {
      if (startsWith(input$fasta,">") == TRUE){
        temp_dir <- tempdir()
		tstamp <- paste(temp_dir,"/input.fasta",sep="")
		cat(input$fasta, file=tstamp)
        input_table <- read.dna(tstamp, format="fasta", as.matrix = FALSE)
        unlink(tstamp)
		input_gene_list <- names(input_table)
		input_datasetname <- input$dataset
		#add datasets from the website dataset option to choose which genes are analyzed.
		#the dataset is assumed to have a folder with the same name as the dataset
		#containing fasta format alignments of each gene in the dataset.
		#Make sure taxa names match across alignments otherwise concatenating will not work.
		#Add more if/else blocks to add more datasets matching your index.html file.
		if(input_datasetname == "testdataset"){
			#use genes listed in test_defined_genes and use the files in the
			#test-dataset folder containing the alignments
			defined_genes <- test_dataset_defined_genes
		} else {
			#should not ever get here, just in case default to test dataset
			defined_genes <- test_dataset_defined_genes
		}

		is_first_gene <- TRUE
 		#foreach gene name
		#   check if it is in the list of our genes
		#   if it is:
		#       change user's sequence name to query_isolate so concatenating works
		#       mafft add to that genes alignment
		#       concatenate with master alignment
		for(i in input_gene_list){
			#check if the user supplied gene name is in our list of allowed genes
			if (i %in% defined_genes){
				gene_aln_file <- paste0(input_datasetname, "/", i, ".aln", "")
				gene_aln <- read.dna(gene_aln_file, format = "fasta")
				input_seq <- input_table[names(input_table)==i]
				names(input_seq) <- c(query_name)
				#add the user's sequence for that gene to the existing alignment generated by mafft
				individual_aln <- mafft(gene_aln, input_seq, "add", method = "localpair", path = mafft_bin, quiet = TRUE)
				if(is_first_gene){
					#Is the first alignment, so use that as the master alignment
					all.al <- individual_aln
					is_first_gene <- FALSE
				} else {
					#otherwise concatenate it with the existing master alignment using phyloch's c.genes()
					temp_al <- c.genes(all.al, individual_aln)
					all.al <- temp_al
				}
			} else {
				#gene name given by the user is not in the dataset, return an error
				return(FALSE)
			}
		}
        return(all.al)
      }else{
        return(FALSE)
      }
    }
    })

# Matching the distances from the user interface (<select id="distance" style="width:300px">)
  dist.genoid <- reactive({
    all.dist <- dist.dna(alin(),input$model)
    return(all.dist)
  })

# Importing data into poppr/adegenet
  data.genoid <- reactive({
    gen <- DNAbin2genind(alin())
    #Adding colors to the tip values according to the clonal lineage
    popnames <- get_last_substring(indNames(gen), "_")
    pop(gen) <- popnames
    gen$other$tipcolor <- pop(gen)
    name_char_end <- nchar(indNames(gen)) - nchar(popnames) - 1
    indNames(gen) <- substr(indNames(gen), 1, name_char_end)
    gen$other$input_data <- indNames(gen)[pop(gen) %in% "query"]
    ngroups   <- length(levels(gen$other$tipcolor))
     ########### MICROBE-ID customization  ############
     # Change these colors to represent the groups defined in your data set.
     #
     defined_groups <- c("blue")
     #
     # Change heat.colors to whatever color palette you want to represent
     # submitted data.
     #
     input_colors   <- rainbow(ngroups - length(defined_groups))
     #
     ##################################
    levels(gen$other$tipcolor) <- c(defined_groups, input_colors)
     gen$other$tipcolor <- as.character(gen$other$tipcolor)
    return(gen)
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

  # Distance tree with bootstrap
  boottree <- reactive({
    if (input$boot > 1000){
      return(1000L)
    } else if (input$boot < 10){
      return(10L)
    }
    set.seed(seed())
      if (input$tree == "upgma"){
        tre <- upgma(dist.genoid())
        bp <- boot.phylo(tre,alin(),function(x) upgma(dist.dna(x,input$model)),B=input$boot)

  	  } else if (input$tree == "njs") {
        tre <- njs(dist.genoid())
	    bp <- boot.phylo(tre,alin(),function(x) njs(dist.dna(x,input$model)),B=input$boot)
  	  } else {
        tre <- nj(dist.genoid())
        bp <- boot.phylo(tre,alin(),function(x) nj(dist.dna(x,input$model)),B=input$boot)
      }
    tre$node.labels <- round(((bp / input$boot)*100))
	#collapse short branches into polytomies
	tre <- di2multi(tre)
	tre$tip.label <- paste(tre$tip.label, as.character(pop(data.genoid())))
    if (input$tree=="nj" || input$tree=="njs"){
      tre <- phangorn::midpoint(ladderize(tre))
    }
    return(tre)
  })

# Minimum spanning network
  msnet <- reactive ({
    msn.plot <-poppr.msn(data.genoid(), dist.genoid(), palette=rainbow, showplot = FALSE)
    return(msn.plot)
  })

############ MICROBE-ID customization ############
# The following lines of code communicate with the user interface to
# plot the outputs from the processes in the server script.

# Validating query as FASTA file
output$validateFasta <- renderText({
  if (is.null(alin())){
    return("")
  } else if (alin() == FALSE){
    return("ERROR: Invalid FASTA file")
  }
  })

## Distance Tree 	(<div id="distPlotTree" class="span6 shiny-plot-output">)
  output$distPlotTree <- renderPlot({
    if (is.null(alin())){
      plot.new()
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) +
      text(x = 0.5, y = 0.9, "No FASTA data has been input.", cex = 1.6, col = "white")
    } else if (alin() == FALSE){
      plot.new()
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) +
        text(x = 0.5, y = 0.9, "Invalid FASTA input.", cex = 1.6, col = "white")
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
    if (is.null(alin())){
      plot.new()
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) +
        text(x = 0.5, y = 0.9, "No FASTA data has been input.", cex = 1.6, col = "white")
    } else if (alin() == FALSE){
      plot.new()
      rect(0, 1, 1, 0.8, col = "indianred2", border = 'transparent' ) +
        text(x = 0.5, y = 0.9, "Invalid FASTA input.", cex = 1.6, col = "white")
    } else {
      plot.minspan(data.genoid(), msnet(), gadj=c(slider()), ind = data.genoid()$other$input_data)
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
      plot.tree(boottree(), type = input$tree, tip.col=as.character(unlist(data.genoid()$other$tipcolor)))
      dev.off()
    })

## Minimum spanning network in PDF format (<a id="downloadPdfMst"  class="btn btn-info shiny-download-link")
  output$downloadPdfMst <- downloadHandler(
    filename = function() { paste0("min_span_net", '.pdf')} ,
    content = function(file) {
      pdf(file, width=11, height=8.5)
      set.seed(seed())
      plot.minspan(data.genoid(), msnet(), gadj=c(slider()), ind = data.genoid()$other$input_data)
      dev.off()
    })
#EOF
})
