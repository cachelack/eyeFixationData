# eyeFixationData

Introduction:

The following is an overview of the R code associated with 
the article "Markov models for ocular fixation locations in 
the presence and absence of colour" by Adam B Kashlak, 
Eoin Devane, Helge Dietert, and Henry Jackson from the
University of Cambridge.  Any questions or comments can 
be directed at Adam B Kashlak (ak852@cam.ac.uk).

Note that most all of the functions in the R code 
begin with the prefix 'rss'.  For the curious, the 
reason is because most of this code was written in 
2015 RSS Statistical Analytics Challenge.
(see https://rsschallenge.wordpress.com/)

Files Included:

The R code: 
  eyeMarkovCode4JRSSC.r
The data table:
  data.csv
The directory of images:
  images/


Data:

To access the data from the R code, the user must set the 
variable 'dataDir' at the top of the .r file.  This currently 
defaults to the current working directory "./".  The data table is 
contained in ./data.csv.  The image files can be 
found in the subdirectories of ./images/.  For the
sake of statistics, the image files are not necessary.  However,
if they are absent, then many of the display functions in the
R code will fail to run correctly.

Data is available on GitHub at 
https://github.com/cachelack/eyeFixationData.git

Code Examples:

The R code is a bit cumbersome as it has not been cleaned into
a nice library.  Here are some examples of functions to run
to reproduce results from the article.

***
* To display photos with plotted fixations with 
  'photo_number' in {1,2,...,60} and 'type' in 
  {'normal','abnormal','grayscale'}

rss_plotImageData2( photo_number,type )

***
* To extract some subset of the data.csv table with
  image     imag in {1,...,60} 
  condition type in {"normal","grayscale","abnormal"}
  subject   subj in {1,...,10}
  fixation  fixt in {1,...,26}

rss_getData <- function( imag=F, type=F, subj=F, fixt=F )

* Similarly, to extract data about the transitions, which 
  is the Euclidean distance between successive fixations.

rss_getSacData <- function( imag=F, type=F, subj=F, fixt=F )

***
* To recreate Figure 6,

rss_plotDurations()

***
* From Section 2 of the article:
  To compute Bayes factors for model selection and choose 
  an optimal clustering, use rss_computeAllBayes( ).
  For example, so analyze the first 2 photos over all 
  three colour schemes, run

rss_computeAllBayes( images=1:2 )

* This function basically calls rss_mmBayesFactor( ) on 
  all of the specified images.  This can be called individually
  as follows

out = rss_mmBayesFactor( 1, 'normal' )

  It produces all of the necessary output for further analysis
  as in the ROC curve from Figure 3



