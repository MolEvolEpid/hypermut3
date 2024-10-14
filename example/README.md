## Example data, output, and plot

This example data comes from sequences from [Patient ID 99722](https://www.hiv.lanl.gov/components/sequence/HIV/search/patient.comp?pat_id=99722&id=1aa26c820e35983d15a94341b42d1a40) from the LANL HIV database.

To run the example data through hypermut, navigate to this directory and run the command:

```
bash run_hypermut.sh
```

Assuming that you have R and tidyverse on your computer, to generate the example plot, you can run the command:

```
Rscript example_plot.R
```

The example plot will be in `example.pdf`. You will see that there is a big difference in the cumulative number of matches for the resolved consensus sequnce, which is the one that conatins multistate characters. 
