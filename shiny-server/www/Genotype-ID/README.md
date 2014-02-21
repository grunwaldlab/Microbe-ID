Genotype-ID
===========

Creation of a Shiny server to develop population genetic tool ocurred to us while working on [poppr] (https://github.com/poppr/poppr) on the winter of 2013.

The goals of the Genotype-ID are
* To create an easy framework for the reconstruction and identification on SSR data
* To use [poppr] (https://github.com/poppr/poppr) in a highly efficient way, in order to create web-apps for biologists/researchers everywhere
* To generate a genotyping webpage for species of *Phytopthora*

Species implemented so far
--------------------------
+ *Phytophthora ramorum*



Installation
----------------------

**Genotype-ID** is based on R, it uses:
* The latest [R](http://cran.r-project.org/) version (Version 2.15.3 is currently supported)
* [poppr] (https://github.com/poppr/poppr) package and its dependencies. 
  ***Note:*** Right now, poppr must be in bleeding edge version for the NJ to work. In this case, refer to [poppr manual](http://grunwaldlab.cgrb.oregonstate.edu/primer-population-genetic-analyses-r/installation) to install it. 
* [Shiny server] (https://github.com/rstudio/shiny-server) and its dependencies

After installing the requirements, be sure that you have a folder named **shiny-server** on ``/var/``. Download the files into ``/var/shiny-server/`` using ``git pull``:

``$ git pull https://github.com/Tabima/Genotype-ID/``

After downloading the files, start the shiny-server

``$ sudo shiny-server &``

App Deployment
-----------------------

And go to your preferred browser (Genotype-ID is known to have issues with Internet Explorer, so we recommend using [Google Chrome](https://www.google.com/intl/es/chrome/browser/?hl=es) or, better yet, [Firefox](http://www.mozilla.org/en-US/firefox/new/#download-fx)) and go to the server port address:

``http://localhost:3838/``

The browser will show you a link to a folder called "Genotype-ID/", if you click on it, the application will deploy on shiny server and voila! there is your web-app.




