#!/usr/bin/perl

use CGI;
#use CGI.pm module

$q=new CGI;

$jobname=$q->param('job_name');
$fasta=$q->param('fasta_file');

#get the parameter from name field
# and store in $value variable.

open (FASTA, ">query.fasta") or die ("Could not open file");
print FASTA $fasta;
close FASTA;
print $q->header;
&html;

sub html{
print <<EOF
<!DOCTYPE html>
<html lang="en"><head>
<meta http-equiv="content-type" content="text/html; charset=UTF-8">
    <meta charset="utf-8">
    <title>Microbe-ID | BLAST Results: $jobname</title>
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="description" content="">
    <meta name="author" content="root" >

    <!-- Le styles -->
    <link href="http://twitter.github.io/bootstrap/assets/css/bootstrap.css" rel="stylesheet">
    <style>
      body {
        padding-top: 60px; /* 60px to make the container go all the way to the bottom of the topbar */
      }
    </style>
    <link href="Bootstrap_files/bootstrap-responsive.css" rel="stylesheet">
        <link href="Bootstrap_files/phytophthora-lab.css" rel="stylesheet">

    <!-- Fav and touch icons -->
    <link rel="apple-touch-icon-precomposed" sizes="144x144" href="http://twitter.github.com/bootstrap/assets/ico/apple-touch-icon-144-precomposed.png">
    <link rel="apple-touch-icon-precomposed" sizes="114x114" href="http://twitter.github.com/bootstrap/assets/ico/apple-touch-icon-114-precomposed.png">
      <link rel="apple-touch-icon-precomposed" sizes="72x72" href="http://twitter.github.com/bootstrap/assets/ico/apple-touch-icon-72-precomposed.png">
                    <link rel="apple-touch-icon-precomposed" href="http://twitter.github.com/bootstrap/assets/ico/apple-touch-icon-57-precomposed.png">
                                   <link rel="shortcut icon" href="http://phytophthora-id.org/img/customIcon.png">
  </head>
 <body>
 <div class="navbar navbar-inverse navbar-fixed-top">
      <div class="navbar-inner">
        <div class="container">
          <button type="button" class="btn btn-navbar" data-toggle="collapse" data-target=".nav-collapse">
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="brand" href="/index.html">Microbe-ID</a>
          <div class="nav-collapse collapse">
            <ul class="nav">
              <li class="/index.html"><a href="/index.html">Home</a></li>
              <li><a href="/seq-id.html">Sequence ID</a></li>
              <li><a href="/geno-id.html">Genotype ID</a></li>
              <li><a href="/about.html">About</a></li>            </ul>
          </div><!--/.nav-collapse -->
        </div>
      </div>
    </div>

    <div class="container">
	<h1>Blast results</h1>
 	<p><strong>Selected jobname is:</strong>$jobname </p>
	<hr></hr>
   	<h2 class="text-info"> Result table </h2>
EOF
;
##########
# IMPORTANT NOTE: Change the path of the blastn executable to where you have it installed in your machine
# Also, don't forget to create the database for your BLAST search using the makeblastdb command that comes with the blast binaries.
# Change the path of said database in the -db argument
###########
$command="/usr/local/bin/ncbi-blast-2.2.28+/bin/blastn -query query.fasta -db /data/www/cgi-bin/phyto_id2/Phytophthora_id/Phytopthora-ID-ITS-version-1.7.fasta -num_alignments 20 -num_descriptions 20 -html";
system("$command");
	
print <<EOF 

    </div> <!-- /container -->

    <!-- Le javascript
    ================================================== -->
    <!-- Placed at the end of the document so the pages load faster -->
<script src="Bootstrap_files/jquery.js"></script>
	<script src="Bootstrap_files/bootstrap-transition.js"></script>
	<script src="Bootstrap_files/bootstrap-alert.js"></script>
	<script src="Bootstrap_files/bootstrap-modal.js"></script>
	<script src="Bootstrap_files/bootstrap-dropdown.js"></script>
	<script src="Bootstrap_files/bootstrap-scrollspy.js"></script>
	<script src="Bootstrap_files/bootstrap-tab.js"></script>
	<script src="Bootstrap_files/bootstrap-tooltip.js"></script>
	<script src="Bootstrap_files/bootstrap-popover.js"></script>
	<script type="text/javascript" src="http://jzaefferer.github.com/jquery-validation/jquery.validate.js"></script>
 </body>
      <footer class="footer">
      <div class="row-fluid">
 		<p> Add some text here <p>
     </div>
      </footer>
      </html>
EOF
;}
;
