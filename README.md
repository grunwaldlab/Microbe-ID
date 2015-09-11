Microbe-ID
===============

This site provides all code and applications necessary for creating a website for species identification or genotyping. A working implementation of Microbe-ID for plant pathogens in the genus *Phytophthora* can be found at [Phytophthora-ID](http://phytophthora-id.org).  

Components
------

To install Microbe-ID on your server you need the following components and applications:

- Access to a linux server to host the site.
- A website based on [bootstrap](http://getbootstrap.com) to host your version of Microbe-ID. See our example [here](./Bootstrap_files).
- A self-hosted [shiny](http://www.rstudio.com/shiny/) R server session running on your server.
- Two [shiny R scripts](./shiny-server/www/Readme.md) for server and user interface (UI), respectively.
- Your custom fasta database for **Sequence-ID**.
- A custom file input for each of the genotyping modules (`Genotype-ID`, `MLST-ID` or `Binary-ID`, each has a description and example in the `shiny-server` folder).
- Further custom modules or pages can be added as needed.

Customization
-------

Follow the following steps to cretae your custom website:

- Customize each html file to your needs including for example ['index.html'](./index.html), ['seq-id.html'](./seq-id.html), ['geno-id.html'](./geno-id.html), ['MLST.html'](./MLST.html), ['binary-id.html'](./binary-id.html), ['about.html'](./about.html), remove unnecesary files and add any additional pages for your site. Instructions for customizations are provided inside each html document as comments starting with '<!--MICROBE-ID customization: ...'.
- Customize each server.R and ui.R file required by shiny for each applitaiton found under [shiny-server/www](./shiny-server/www). 
- Do not modify any other files including ['validate.js'](./valdiate.js).
