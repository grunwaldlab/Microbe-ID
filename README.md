Microbe-ID
===============

Microbe-ID is a toolbox that can readily be customized for any taxonomic group or species for sequence based identification of species or marker-based identification of strains. Microbe-ID aims to:

1. Outline a method for development of a species and strain identification website using existing open source, online resources

2. Integrate existing population genetic analyses and visualizations into the genotype identification tool in order to improve genotype diagnostic capabilities.

This site provides all code and applications necessary for creating a website for species identification or strain genotyping.

## Components

Microbe-ID is comprised by two main applications: **Sequence-ID** and **Genotype-ID**. **Sequence-ID** allows identification of species by blasting a query sequence for any locus of interest against a custom reference sequence database using BLAST.  **Genotype-ID**  allows placement of an unknown genotype in either a minimum spanning network or dendrogram with bootstrap support comparing the queue against an user-created reference database using R and a suite of R packages designed for molecular markers. **Microbe-ID** groups these two applications in a HTML5 framework based on [Bootstrap](http://getbootstrap.com) that links both applications into a simple, web-ready site to be deployed locally or in a server.

### Files included in this repository

#### Front end HTML files
- [`index.html`](./index.html): Front page of **Microbe-ID**.
- [`seq-id.html`](./seq-id.html): Front page of **Sequence-ID**. Provides the user with forms to submit the query sequences to be processed by BLAST to identify similarity in a custom-well curated database. This page provides the informatiom and connections to the `blast.cgi` script, which runs BLAST using PERL.
- [`geno-id.html`](./geno-id.html): Front page of **Genotype-ID**. Provides the user with examples, forms and instructions in how to submit **SSR/Microsatellite** queries to be processed by R using shiny. This page provides a iFrame to [`shiny-server/www/Genotype-ID/www/index.html`](./shiny-server/www/Genotype-ID/www/index.html).
- [`MLST.html`](./MLST.html): Front page of **MLST-ID**. Provides the user with examples, forms and instructions in how to submit **Multi Locus Seqeunce Type (MLST)**  queries to be processed by R using shiny. This page provides a iFrame to [`shiny-server/www/MLST/www/index.html`](./shiny-server/www/MLST/www/index.html).
- [`binary-id.html`](./binary-id.html): Front page of **Genotype-ID**. Provides the user with examples, forms and instructions in how to submit **Codominant/Binary** queries to be processed by R using shiny. This page provides a iFrame to [`shiny-server/www/Binary-ID/www/index.html`](./shiny-server/www/Binary-ID/www/index.html).
- [`about.html`](./about.html): Contact page.

#### Front end HTML folders
- [`Bootstrap_files/`](./Bootstrap_files/): Folder that contains all the files necessary for the [Bootstrap front-end framework](getbootstrap.com/) in which **Microbe-ID** was designed and built from.
- [`cgi-bin/`](./cgi-bin/): Folder that contains the `blast.cgi` script for **Sequence-ID**. `blast.cgi` is a PERL-CGI script that runs `BLAST` in a custom database of curated sequences.
- [`files/`](./files/): Folder in which additional files can be added.
- [`img/`](./img/): Folder that contains images used in **Microbe-ID** webpages.

#### Shiny applications
- [`shiny-server/www/`](./shiny-server/www/): Contains all the folders and files to be used by the shiny server.
  - [`Binary-ID/`](./shiny-server/www/Binary-ID/): Contains all the folders and files to be used by the shiny server for the deployment of **Binary-ID**.
    - [`Aeut.txt`](./shiny-server/www/Binary-ID/Aeut.txt): Curated database of AFLP molecular markers used as references in **Binary-ID**. This file serves as input for `server.R`.
    - [`server.R`](./shiny-server/www/Binary-ID/server.R): File that contains all the commands required by R to execute **Binary-ID** by R. This file is handled by the shiny server.  More information available on the **Shiny Arquitecture** region.
    - [`www/`](./shiny-server/www/Binary-ID/www): Contains the files necesary to deploy the **user interphase** used by **Binary-ID**.
      - [`AFLP_Aphanomyces.xlsx`](./shiny-server/www/Binary-ID/www/AFLP_Aphanomyces.xlsx): The `AFLP_Aphanomyces.xlsx` file contains a template to format the queries as an example for the end users of **Binary-ID**.
      - [`index.html`](./shiny-server/www/Binary-ID/www/index.html): User interphase file for **Binary-ID**. More information available on the **Shiny Arquitecture** region.
  - [`Genotype-ID/`](./shiny-server/www/Genotype-ID/): Contains all the folders and files to be used by the shiny server for the deployment of **Genotype-ID**.
    - [`Ramorum_ssr.csv`](./shiny-server/www/Genotype-ID/Ramorum_ssr.csv): Curated database of SSR/Microsatellite molecular markers used as references in **Genotype-ID** for *Phytopthora ramorum*. This file serves as input for `server.R`.
    - [`server.R`](./shiny-server/www/Genotype-ID/server.R): File that contains all the commands required by R to execute **Genotype-ID** by R. This file is handled by the shiny server. More information available on the **Shiny Arquitecture** region.
    - [`www/`](./shiny-server/www/Genotype-ID/www): Contains the files necesary to deploy the **user interphase** used by **Genotype-ID**.
      - [`SSR_Example_Data_ramorum.xlsx`](./shiny-server/www/Genotype-ID/www/SSR_Example_Data_ramorum.xlsx): The `SSR_Example_Data_ramorum.xlsx` file contains a template to format the queries as an example for the end users of **Genotype-ID**.
      - [`index.html`](./shiny-server/www/Genotype-ID/www/index.html): User interphase file for **Genotype-ID**. More information available on the **Shiny Arquitecture** region.
  - [`MLST/`](./shiny-server/www/MLST/): Contains all the folders and files to be used by the shiny server for the deployment of **MLST-ID**.
    - [`server.R`](./shiny-server/www/MLST/server.R): File that contains all the commands required by R to execute **MLST-ID** by R. This file is handled by the shiny server. More information available on the **Shiny Arquitecture** region.
    - [`www/`](./shiny-server/www/MLST/www): Contains the files necessary to deploy the **user interphase** used by **MLST-ID**.
      - [`example.html`](./shiny-server/www/MLST/www/example.html): HTML page with examples for **MLST-ID** query formatting.
      - [`index.html`](./shiny-server/www/MLST/www/index.html): User interphase file for **MLST-ID**. More information available on the **Shiny Arquitecture** region.
      - [`help.html`](./shiny-server/www/MLST/www/help.html): Help and tutorials for **MLST-ID** use.
    - [`test-dataset/`](./shiny-server/www/MLST/test-dataset): Contains the `FASTA` files used in several tests by **MLST-ID**.

## Requirements

To install Microbe-ID on your server you need the following components and applications:

### General requirements
- Access to a linux server to host the site.
- A website based on [bootstrap](http://getbootstrap.com) to host your version of Microbe-ID. See our example [here](www.microbe-id.org).

### For the R-based apps of Genotype-ID
- A self-hosted [shiny](http://www.rstudio.com/shiny/) R server session running on your server.
- Two [shiny R scripts](./shiny-server/www/Readme.md) for server (`server.R`) and user interface (`ui.R` or `index.html`), respectively.
- A custom file input for each of the genotyping modules (`Genotype-ID`, `MLST-ID` or `Binary-ID`, each has a description and example in the `shiny-server` folder).

### For Sequence-ID
- Your custom FASTA curated database for **Sequence-ID**. This database
- Access to cgi-bin perl applications.
- `BLAST+` suite
- Custom `FASTA` database of curated sequences for the group in interest

## Shiny app architecture

For the R applications deploy successfully, shiny requires two types of files, the **user interphase file** and the **shiny script**.

- **User interphase (ui.R or index.html):** The user interphase file (UI) permits the communication between the user and `R` using forms and reactive functions. TheUI can be defined as the "web-document" and its where the application will prompt the results and receive the input of the user. Examples of UI are available in any of the shiny-app folders of Microbe-ID (look for the `index.html` file)

- **Server script (server.R):** The server script runs the application in R. The server script contains all the information that R will interpret to run the application (code, internal datasets and interpretations of the input/outputs from and to the user interphase file). Examples of server scripts are available in any of the shiny-app folders of Microbe-ID (look for the `server.R` file)

Both `server.R` and `ui.R` files should be in the same directory. If you prefer to use the `index.html` file, create a folder called `www/` and move the `index.html` file to this folder.

## Installation of Genotype ID

### Prerrequisites

- [R](https://www.r-project.org)
- [Shiny Server](http://shiny.rstudio.com/articles/shiny-server.html)
- Shiny app (`server.R` and `ui.R or index.html`)

1. Install `Shiny Server` and `R`. We recommend to follow the instructions carefully from the developers sites.
For more info in `Shiny Server` please read:
  - [Shiny Server professional guide](http://rstudio.github.io/shiny-server/latest/)
  - [Shiny Server installation](http://rstudio.github.io/shiny-server/latest/#installation)

2. Install any additional software.
  - Install additional packages in `R`. We highly recommend to install the packages as a superuser, to make them available to all uses. For any of the applications on **Genotype-ID**, you will need the packages `shiny`,`poppr`, `pegas`, `igraph`, `phangorn`, `gdata`, `XML`, `phyloch` from CRAN.
  - If you are going to use `MLST-ID`, install `mafft`. `mafft` is a multiple sequence alignment program using fourier algorithms. `MLST-ID` uses `mafft` to align and creates multiple alignments for each of the MLST.
  - If you are planning on using **Sequence-ID**, install `BLAST`. `BLAST` is used by **Sequence-ID** to identify sequence data by similarity to a custom well-curated database.  For more information in `BLAST`, go to the [main page](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs) or to [downloads and installation](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

3. Clone this repository in a folder accessible by Apache (`cgi-bin` integration is required to use Sequence-ID)

4. Move the contents of `shiny-server/www` folder app to the Shiny Server directory on your local server (`/shiny-server/`) and restart your local `shiny-server`. If the installation is correct, the modules will be visible and executable at `127.0.0.1:3838`.

5. If you are installing **Sequence-ID**:

  - Move the `blast.cgi` file from the repository `cgi-bin/` folder into your apache `cgi-bin` folder.
  - Create the database for your custom `BLAST` search using the `makeblastdb` command that comes with the blast binaries.
  - Change the path of the `blastn` executable and the `BLAST` database in the `blast.cgi` file.

6. Integrate the **Genotype-ID** modules into your custom webpage by changing the paths on the `frames` in the ['geno-id.html'](./geno-id.html), ['MLST.html'](./MLST.html), ['binary-id.html'](./binary-id.html) files.


## Customization

### Genotype-ID

- Customize each html file to your needs including for example ['index.html'](./index.html), ['seq-id.html'](./seq-id.html), ['geno-id.html'](./geno-id.html), ['MLST.html'](./MLST.html), ['binary-id.html'](./binary-id.html), ['about.html'](./about.html). Instructions for customizations are provided inside each html document as comments starting with `<!--MICROBE-ID customization: ...>`. Please follow the directions of these tags carefully.
- Remove unnecessary files and add any additional pages for your site.
- Images and additional files can be added to the `img/` and `files/` folders respectively. Images or files can be added to any of the `html` pages using HTML5. More information on HTML5 image tags [here](http://www.w3schools.com/html/html_images.asp)
- Customize each `server.R` and `ui.R/index.html` file required by shiny for each application under [shiny-server/www](./shiny-server/www).
- Do not modify ['validate.js'](./validate.js) or the `Bootstrap_files/` folder as they contain fundamental scripts for the website.
- In case the user wants to modify the `.gitignore` file, please do not remove the files present originally.

### Sequence-ID

- Move the `blast.cgi` file from the repository `cgi-bin/` folder into your apache `cgi-bin` folder.
- Create the database for your custom `BLAST` search using the `makeblastdb` command that comes with the `blast` binaries.
- Change the path of the `blastn` executable and the `BLAST` database in the `blast.cgi` file.

# Example implementation

A working implementation of Microbe-ID for plant pathogens in the genus *Phytophthora* can be found at [Phytophthora-ID](http://phytophthora-id.org). See the sample code on [github](https://github.com/grunwaldlab/phytophthora_id)  for this implementation.
