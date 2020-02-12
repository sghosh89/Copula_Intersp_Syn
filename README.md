# Repository of analyses: A new approach to interspecific synchrony in population ecology using tail association

Shyamolina Ghosh, University of Kansas 

Lawrence W. Sheppard, University of Kansas

Philip C. Reid, University of Plymouth

Daniel C. Reuman, University of Kansas

## Introduction
This repository records the complete workflow that produced the paper from the data. All 
analyses can be reproduced and the paper and supporting information recompiled (see below).

## How to compile
Knit makefile.Rmd using R markdown. If all dependencies are in place (see next section) 
this should re-compute all analyses from data to paper, resulting in three pdfs: 
**MainText.pdf** (the main text of the paper), **SuppMat.pdf** (the 
supporting information file for the paper), and **makefile.pdf** (notes on the 
compilation process - can be useful for error mitigation in the event of failure).

The knit may take quite some time. Subsequent knits, if any, 
can be faster because packages will be installed (see below) and because intermediate 
results are cached. If you try to knit MainText.Rmd or SuppMat.Rmd directly, you may have some 
success, but cross-document references and other features will fail so this is not recommended.
To compile the documents from the command line, use the following: Rscript -e "library(knitr); knit('makefile.Rmd')".

## Dependencies

### Core software dependencies
   - R 
   - R markdown
   - R studio
   - latex 
   - bibtex

### Data dependencies
Datasets are not included in the Data folder, and need to be obtained and put there for the code
to run. A dataset which includes the plankton data we used as a subset can be obtained from the 
Dryad Digital Repository https://doi.org/10.5061/dryad.rq3jc84. We do not have the rights
to release the aphid data, so those data are not in the repository. 
The aphid data came from the 
Rothamsted Insect Survey (RIS) of Rothamsted Research (https://www.rothamsted.ac.uk/insect-survey). The plankton data came from the Continuous Plankton
Recorder (CPR) dataset of the Marine Biological Association of the UK (https://www.cprsurvey.org/). Both organizations 
have clear policies for sharing data on their websites. James Bell (james.bell@rothamsted.ac.uk) 
is our contact at RIS and P. Chris Reid (pchrisreid@googlemail.com) is our contact at CPR. 
If a user obtains written permission from these organizations then we will be happy to provide 
these datasets in the format expected by repository code.

### Dependencies on the R checkpoint package

Code uses the R *checkpoint* package. This is set up in the master file makefile.Rmd in the 
R chunk *checkpoint_chunk*, which contains the following line of code specifying a date :

checkpoint("2019-01-01",checkpointLocation = "./")

The checkpoint package then automatically scans through other files looking for other required R 
packages. It then downloads and installs the newest versions of those packages available on the 
given date. This helps ensure that re-compiling the document uses exactly the same code that was 
originally used, in spite of package updates and other changes. This can take some time on first 
run (you are warned) but it is faster on subsequent runs because the packages are already 
installed. This also means that R package dependencies should only be the checkpoint package, 
since that package should scan for other packages and install them locally. Quite a few MB disk 
space are used.

### Dependencies on pandoc
The open source program pandoc converts documents from one format to another. 
Here, the knitr package uses it to convert the markdown files into latex format so that 
they can then be turned into PDF files. Installers for multiple operating systems are available 
here: https://pandoc.org/installing.html.

### Dependencies on pdflatex
The makefile makes a system call to pdflatex, so software supporting that needs to be installed:
 - On Windows, you can use Miktex (https://miktex.org/howto/install-miktex),
 - On Linux, install latex (e.g., sudo apt-get install texlive), and
 - On Mac, use the MacTeX installer (http://www.tug.org/mactex/)

### Additional dependencies?
If you find additional dependencies were needed on your system, please let us know: 
reuman@ku.edu. The compilation process was tested by Ghosh on Ubuntu 16.04 and by 
Reuman on a similar computing setup. It has not been 
tested on Mac. We have endeavored to list all dependencies we can think of above, but we have 
only compiled on our own machines, so we cannot guarantee that additional dependencies will not 
also be needed on other machines, even after data are included (see above). This repository is 
intended to record a workflow, and is not designed or tested for distribution and wide use on 
multiple platforms. It is not guaranteed to work on the first try without any hand-holding on 
arbitrary computing setups.

## Intermediate files:
Knitting the makefile automatically produces a lot of 'intermediate' files. Files ending in .tex 
are the converted documents from .Rmd including all the R code output and the rest (files ending 
.log, .aux, .lof, .lot, .toc and .out ) are intermediate files that pdflatex uses to keep track 
of various parts of the document. Some of these can be useful for diagnosing problems, if any.

## Acknowlegements:
We thank the many contributors to the large datasets we used; D. Stevens and P. Verrier for 
data extraction; and Joel E. Cohen, Lauren Hallet, and Jonathan Walter for helpful 
suggestions. We thank James Bell of the Rothamsted Insect Survey (RIS). The RIS, a UK Capability,
is funded by the Biotechnology and Biological Sciences Research Council under the Core
Capability Grant BBS/E/C/000J0200. SG, LWS and DCR were partly funded by US National Science
Foundation grants 1714195 and 1442595 and the James S McDonnell Foundation. 
Any opinions, findings, and conclusions or recommendations expressed in this 
material are those of the author(s) and do not necessarily reflect the views of the National 
Science Foundation or the McDonnell Foundation.











