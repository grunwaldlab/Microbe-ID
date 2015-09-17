# Microbe-ID shiny R modules

We created several `shiny R` modules for identification of microbial taxa using different molecular markers such as SSR/microsatellites, MLST, and dominant or binary markers (AFLP, RFLP, etc.).

Each of the modules are open source, fully customizable and easy to use, in order to provide the researcher with a modular web-application that permits flexible strain identification for virtually any microbial species.

## General Components

Quoting the shiny Rstudio tutorial:

```
Shiny apps have two components:

* a user-interface script: The user-interface (ui) script controls the layout and appearance of your app.

* a server script: The script contains the instructions that your computer needs to build your app.

In our cases, for each of the the components we have a server script (`server.R`) and an HTML based user interphase (`www/index.html`)

The user interphase communicates with the R application by using forms. Each of the forms will contain the query information that will be passed to the R app and be analysed with the reference datasets. Each of the components is fully customizable and contain instructions to customize them efficiently

```

## Installation


**Genotype-ID** is based on R, it uses:
* The latest [R](http://cran.r-project.org/) version
* [poppr] (https://github.com/poppr/poppr) package and its dependencies.
* [Shiny server] (https://github.com/rstudio/shiny-server) and its dependencies

After installing the requirements, be sure that you have a folder named **shiny-server** on ``/var/``. Download the files into ``/var/shiny-server/``.

After downloading the files, start the shiny-server

``$ sudo shiny-server &``

App Deployment
-----------------------

And go to your preferred browser (Genotype-ID is known to have issues with Internet Explorer, so we recommend using [Google Chrome](https://www.google.com/intl/es/chrome/browser/?hl=es) or, better yet, [Firefox](http://www.mozilla.org/en-US/firefox/new/#download-fx)) and go to the server port address:

``http://localhost:3838/``

The browser will show you a link to a folder called "Genotype-ID/", if you click on it, the application will deploy on shiny server and voila! there is your web-app.
## Genotype-ID

Genotype-ID is a strain identification tool based on SSR. Genotype-ID uses `poppr` to construct a **Bruvo's distance dendrogram** (using UPGMA/NJ) or a **Minimum spanning network** of an SSR query against a well-curated database of reference SSR's.

### Requirements

* `poppr` package

### Database:

Tab-delimited text file. Column names are loci names, rownames are samples. Each of the alleles is separated by a */*. The second column contains the reference population/lineage/code to group reference strains into discrete groups.

Here's an example of the `reduced_database.txt.csv` file that is used by *Genotype-ID* as reference:

|Sample|Lineage|PrMS6|PrMS9c3|ILVOPrMS132/PrMS39a|PRMS45|ILVOPrMS133/PrMS43ab|KI18|KI64|ILVOPrMS131|
|---|---|---|---|---|---|---|---|---|---|
|PR-12-044|NA1|165/168|216/226|216/226|372/489|130/246|219/275|342/379|150/234|
|PR-11-015|NA1|165/168|216/226|216/226|372/481|130/250|219/275|342/379|154/222|
|PR-11-008|NA1|165/168|216/226|216/226|368/485|130/250|219/275|342/379|150/226|

### Input:

At the user interphase, a text form will be filled containing the same loci for different queries:

```
sample1	???	169/172	221/231	136/252	172/192	373/489	224/278	347/383	150/234
sample2	???	169/172	221/221	152/158	0/0	173/173	226/226	349/365	278/342
```
### Results:

* Customizable minimum spanning network
* UPGMA/NJ dendrogram

## MLST-ID

MLST-ID is a MultiLocus Sequence Type Id app created in shiny R. It uses the information of multiple `FASTA` files of reference organisms to provide clear and fast identification of microbial species using MLST.

### Requirements:

*R packages*:

* `ape`
* `phangorn`
* `XML`
* `phyloch`

*Additional software*

* [MAFFT](http://mafft.cbrc.jp/alignment/software/)

### Database:

The MLST-ID uses a `FASTA` reference per gene/locus of interest:

`Gene_A.fasta`:

```
>0310_WAS
CGCGGCAGTCGCCGTCGCCGAGGCGCCGGCCAAGGCATACAACCCTCTCTTCATCTACGGCGACTCGGGCCTCGGCAAGACCCACCTGCTGCACGCCATC
>0202_ECU
CGCGGCAGTCGCCGTCGCCGAGGCGCCGGCCAAGGCATACAACCCTCTCTTCATCTACGGCGACTCGGGCCTCGGCAAGACCCACCTGCTGCACGCCATC
```
### Input:

The input is a text form that requires a `FASTA` input. The database which the query will be compared needs to be specified in the FASTA identifier (>Gene_A):

```
>Gene_A
CGCGGCAGTCGCCGTCGCCGAGGCGCCGGCCAAGGCATACAACCCTCTCTTCATCTACGGCGACTCGGGCCTCGGCAAGACCCACCTGCTGCACGCCATC
```
The query will be aligned with the database of the reference sequences for the gene of interest using `MAFFT`, and then the selected genetic distance will be calculated to construct the dendogram/minimum spanning network of the query/database.

### Results:

* Customizable minimum spanning network
* UPGMA/NJ dendrogram

### Binary-ID

Binary-ID uses codominant molecular markers (coded as 1 and 0) and different genetic distances to reconstruct a distance dendogram with bootstrap support values and a minimum spanning network.

### Requirements:

* `poppr` R package

### Database:

A dataset of AFLP loci of *Aphanomyces eumithes* of two different populations, as included in `poppr` was used a an example dataset. `Aeut.txt`: Column names are each of the AFLP site profiles. Row names are sample names. Second column are the populations of the reference samples:

|Ind|Pop|1|2|3|4|5|6|7|
|---|---|---|---|---|---|---|---|---|
|1|Athena|1|1|1|0|0|0|0|
|2|Athena|0|1|1|0|0|0|0|
|3|Mt. Vernon|1|1|0|0|0|0|0|
|4|Mt. Vernon|1|1|0|0|0|0|0|
|5|Mt. Vernon|1|1|0|0|0|0|0|

### Input:

At the user interphase, a text form will be filled containing the same loci for different queries:

```
Ind_1	Queue	1	0	1	0	1	0	0
Ind_2	Queue	1	0	1	0	0	0	0
```
### Results:

* Customizable minimum spanning network
* UPGMA/NJ dendrogram
