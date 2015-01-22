Nystrom Bestiary v2.0  
=====================
written by [Alex Gittens](http://thousandfold.net/cz)  
licensed under the [Creative Commons ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/)

It is a collection of code for experimenting with various SPSD Sketches, including Nystrom extensions based on column sampling, Nystrom extensions based on random mixtures of columns, and 'pinched' and 'prolonged' eigensketches (see the review paper ["Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions"](http://arxiv.org/abs/0909.4061) equations (5.11) and (5.12) for the precise definition of these two sketches).

It was used to produce the figures in the paper ["Revisiting the Nystrom Method for Improved Large-scale Machine Learning"](http://arxiv.org/abs/1303.1849) (arXiv preprint link) by Alex Gittens and [Michael Mahoney](http://cs.stanford.edu/people/mmahoney/). In particular, the experimental setup to generate exactly those figures is included.

Send comments, suggestions, and complaints to agittens AT ebay FULLSTOP com

Contents
--------
- the `extensions/` directory contains the implementations of various Nystrom extensions
- the `io/` directory contains the code used to create, load, and process the datasets
- the `datasets/` directory contains the datasets used in the experiments
- the `experiments/` directory contains a set of m-files that actually runs the Nystrom
 extensions on various datasets and stores statistics on the errors and timing
- the `outputs/` directory is used to store the output of the experiments
- the `plots/` directory stores the plots of the timings and errors
- the `auxiliary/` directory contains code needed in computing the extensions
- the `visualization/` directory contains the code used to produce the plots of the
 various timings and errors
- the `misc/` directory contains miscellany (so far, the code to generate the data
for Table 2 in the paper)

Usage instructions
-------
_ALL m files should be run from the base folder, otherwise you'll run into path issues_    
To produce the figures in the paper:  
#####Short story
Ensure that you are in the base directory, `NystromBestiary`, and run the following commands from the Matlab prompt:

    addpath 'auxiliary' 'datasets' 'experiments' 'extensions' 'io' 'misc' 'outputs' 'plots' 'visualization'
    create_bestiary_datasets
    maxNumCompThreads = 1; # if you want accurate timing info
    runall
    visualizeall

#####Long story
1.  add all the subdirectories in this folder to your path
2.  run `create_bestiary_datasets` to generate
 some required distance matrices; this step generates about 1.5Gb of data
3. If you want to have email notifications at the start and end of each
experiment, modify `runall.m` to set the `sendEmails` flag to true and set the 
email-related variables appropriately, then run
`setpref('Internet', 'SMTP_Password', 'youremailpassword')`
at the Matlab command line
4.  run `runall` (or pick individual experiments) in the experiments directory;
 this step generates about 2.7Gb of data
5.  wait several days for the experiments to stop running!
6.  run `visualizeall`

The pdfs will be located in the output directory

See the individual m-files for more details. Make appropriate
modifications to substitute your own datasets.

Attributions:
------------
* jdqrpcg.m is due to [Yvan Notay](http://homepages.ulb.ac.be/~ynotay/) (see the m file for full attribution)  
* notifier.m is due to [Benjamin Krause](http://www.mathworks.com/matlabcentral/fileexchange/authors/74912) (see the m file for full attribution)  

for dataset provenances, see Table 3 in the above mentioned paper (datasets: Abalone, Wine, Spam, Kin8nm, Dexter, Gisette, Enron, Protein, SNPs, HEP, GR, Gnutella)  
two additional datasets, Cranfield and Medline, are from the [Text to Matrix Generator](http://scgroup20.ceid.upatras.gr:8000/tmg/) Matlab Toolbox's website.
