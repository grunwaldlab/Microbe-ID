Microbe-ID
===============

This site provides all code and applications necessary for creating a website for species identification or genotyping. A working implementation of **Microbe-ID** for plant pathogens in the genus *Phytophthora* can be found at [Phytophthora-ID](http://phytophthora-id.org).  

Components
------

To install Microbe-ID on your server you need the following components and applications:

- Your own linux server
- A website based on [bootstrap](http://getbootstrap.com) to host your version of Microbe-ID
- A self-hosted [shiny](http://www.rstudio.com/shiny/) R server session running on your server
- Two shiny R scripts for server and user interface (UI): (https://github.com/grunwaldlab/Microbe-ID/blob/master/shiny-server/www/Readme.md)
- Your custom fasta database for **Sequence-ID**
- A custom file input for each of the genotyping modules (`Genotype-ID`, `MLST-ID` or `Binary-ID`, each has a description and example in the `shiny-server` folder)
- Further custom modules or pages can be added as needed

This site was built solely with open source components.

Designed and built by members of the [Grunwald lab](http://grunwaldlab.cgrb.oregonstate.edu). For questions, please contact [Javier Tabima](caifaz01@gmail.com).
