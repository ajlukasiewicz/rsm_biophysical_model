## RsmA biophysical model 
Code to accompany "Thermodynamic Modeling of RsmA - mRNA Interactions Capture Novel Direct Binding Across the *Pseudomonas aeruginosa* Transcriptome"
Preprint for this work can be found on [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.08.01.606018v1)

## Organization 
Scripts for executing and analyzing model data are provided in [scripts](https://github.com/ajlukasiewicz/rsm_biophysical_model/tree/main/scripts) 

[data](https://github.com/ajlukasiewicz/rsm_biophysical_model/tree/main/data) contains an example of the input format accepted by the model, and two output example files with `binding_sites` and `translation_rate` prefixes

## Usage example
Given an input file in .xls format, `biophysical_model_argparse.py` accepts the following arguments:

```
-i  input file in .xls format
-o  name of output file
-pwm  energy matrix for scoring binding motifs. Accepts CsrA or RsmA_ddG as energy options
-translate  y/n argument for calculating translation rates
-export  file format for exporting model data. Can be "csv" or "excel"
```
An example for running the scripts locally

```
python3 biophysical_model_argparse.py -i ../data/example.xls -o example -pwm RsmA_ddG -translate n -export csv`
```

For running the modeling scripts on compute clusters, the `multiprocess` module is used. 
For these applications you can use `biophysical_model_argparse_MP.py` which uses the same arguments defined above
