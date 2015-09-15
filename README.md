Microbe-ID
===============

Microbe-ID is a toolbox that can readily be customized for any taxonomic group or species for sequence based identification of species or marker-based identification of strains. Microbe-ID aims to:

1. Outline a method for development of a species and strain identification website using existing open source, online resources

2. Integrate existing population genetic analyses and visualizations into the genotype identification tool in order to improve genotype diagnostic capabilities.

This site provides all code and applications necessary for creating a website for species identification or strain genotyping.

## Components

To install Microbe-ID on your server you need the following components and applications:

### General requirements
- Access to a linux server to host the site.
- A website based on [bootstrap](http://getbootstrap.com) to host your version of Microbe-ID. See our example [here](./Bootstrap_files).

### For the R-based apps of Genotype-ID
- A self-hosted [shiny](http://www.rstudio.com/shiny/) R server session running on your server.
- Two [shiny R scripts](./shiny-server/www/Readme.md) for server (`server.R`) and user interface (`UI.r` or `index.html`), respectively.
- A custom file input for each of the genotyping modules (`Genotype-ID`, `MLST-ID` or `Binary-ID`, each has a description and example in the `shiny-server` folder).

### For Sequence-ID
- Your custom fasta database for **Sequence-ID**.
- Access to cgi-bin perl applications.
- `BLAST+` suite
- Custom `FASTA` database of curated sequences for the group in interest

## Shiny app architecture

For the R applications deploy successfully, shiny requires two types of files, the user interphase file and the shiny server.

- User interphase (ui.R or index.html): The user interphase file permits the communication between the user and R using forms and reactive functions. The	UI can be defined as the "web-document" and its where the application will prompt the results and receive the input of the user. Examples of UI are available in any of the shiny-app folders of Microbe-ID (look for the `index.html` file)

- Server file: The server script runs the application in R. The server script contains all the information that R will interpret to run the application (code, internal datasets and interpretations of the input/outputs from and to the user interphase file). Examples of server scripts are available in any of the shiny-app folders of Microbe-ID (look for the `server.R` file)

Both `server.R` and `iu.R` files should be in the same directory. If you prefer to use the `index.html` file, create a folder called `www/` and move the `index.html` file to this folder.

## Installation of Genotype ID

### Prerrequisites

- [R](https://www.r-project.org)
- [Shiny Server](http://shiny.rstudio.com/articles/shiny-server.html)
- Shiny app (`server.R` and `ui.R or index.html`)

1. Install `Shiny Server` and `R`. We recommend to follow the instructions carefully from the developers sites.
For more info in `Shiny Server` please read:
- [Shiny Server professional guide](http://rstudio.github.io/shiny-server/latest/)
- [Shiny Server installation] (http://rstudio.github.io/shiny-server/latest/#installation)

2. Install any additional software.
- Install additional packages in `R`. We highly recommend to install the packages as a superuser, to make them available to all uses. 

If you want to use any of the Genotype-ID modules, you will need the packages `shiny`,`poppr`, `pegas`, `igraph`, `phangorn`, `gdata`, `XML`, `phyloch`
- If you are going to use `MLST-ID`, install `mafft`. `mafft` is a multiple sequence alignment program using fourier algorithms. `MLST-ID` uses `mafft` to align and create a dataset of aligned sites in each of the MLST loci.
- If you are planning on using **Sequence-ID**, install `BLAST`. `BLAST` is used by **Sequence-ID** to identify sequence data by similarity to a custom well-curated database.  For more information in `BLAST`, go to the [main page](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs) or to [downloads and installation](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

3. Clone this repository in a folder accesible by Apache (`cgi-bin` integration is required to use Sequence-ID)

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
- Images and additional diles can be added to the `img/` and `files/` respectively. Images or file can be added to any of the `html` pages using HTML5. More information on HTML5 image tags (here)[http://www.w3schools.com/html/html_images.asp]
- Customize each server.R and ui.R file required by shiny for each application found under [shiny-server/www](./shiny-server/www).
- Do not modify any other files including ['validate.js'](./valdiate.js) and the `Bootstrap_files/` folder.
- In case the user wants to modify the `.gitignore` file, please do not remove the files present originally.

# Example implementation

A working implementation of Microbe-ID for plant pathogens in the genus *Phytophthora* can be found at [Phytophthora-ID](http://phytophthora-id.org). See the sample code on [github](https://github.com/grunwaldlab/phytophthora_id)  for this implementation.
