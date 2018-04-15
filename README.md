# TreeSAP
## A bioinformatics tool for predicting spliced antigenic peptides

Requires Python 3, h2o 3.18.0.2 and pandas.

To install the correct version of h2o with pip enter the following on the command line:

```pip install h2o==3.18.0.2```


## Using TreeSAP

From the command line:

```python3 TreeSAP.py [input] [output]```

Input can be either:
* The name of a file in the TreeSAP directory containing a protein sequence
* A protein sequence in the form of a string

Results are output as a csv file to a file named by the output argument.
