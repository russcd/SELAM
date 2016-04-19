##  A MAKEFILE FOR SELAM. TO USE, RUN THE COMMAND "make" VIA COMMAND LINE ##

TCFLAGS = -ltcmalloc 
GSL_LIBS = $(shell gsl-config --libs)
CFLAGS = $(shell gsl-config --cflags) 

all:
## IF TCMALLOC IS NOT INSTALLED, COMMENT OUT THE FOLLOWING LINE 
## BY PLACING A # IN FRONT OF THE FOLLOWING LINE AND REMOVE '#' ON THE NEXT LINE
	$(LINK.cc) -std=c++0x -O3 SELAM.cpp $(TCFLAGS) $(GSL_LIBS) $(CFLAGS) -o SELAM
#	$(LINK.cc) -std=c++0x -O3 SELAM.cpp $(GSL_LIBS) $(CFLAGS) -o SELAM
	$(LINK.cc) -std=c++0x SELAM_STATS.cpp -o SELAM_STATS


## In many cases, tcmalloc will substantially decrease runtime of SELAM
## we therefore strongly recommend using this library 
## on mac osx it can be installed using homebrew: brew install google-perftools
## on ubuntu distributions you may try: sudo apt-get install google-perftools libgoogle-perftools-dev
