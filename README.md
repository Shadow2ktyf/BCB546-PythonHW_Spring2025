# Python Assignment 
# Please run `PythonHomework_Zheyuan.ipynb`

### Due: 7 May 2025

## Summary

### Document written code and fill in the gaps

Your colleague, Dr. X, has given you a partially written Python script to work with your shared data files. 
You plan to translate cytochrome-b sequences to amino acids for each of 12 species of penguins.
You will then compute some simple summaries of the amino-acid molecular weight and GC content of these sequences and add these values to a DataFrame, which already contains the average adult mass of each species. 
Once you have these values entered in the DataFrame, you will create some plots for data visualization. 

In this script you will find one function and additional lines of code that are working, but not adequately documented. Since Dr. X's experience with Python is different than your own, you may find that some of the suggested functions and approaches are unfamiliar to you; however, you have agreed to thoroughly document your script in order to better understand your work.

Additionally, Dr. X has written some comments and pseudocode that outline other components of the script that you agreed to complete. Since you are a conscientious scientist, you will also carefully document your own code.

### To complete the assignment and analysis, create a single Jupyter notebook that documents the entire workflow. 

## Details

* Create a **_new_** repository for this assignment and name it `BCB546-PythonHW_Spring2025`. The repository you submit should _only_ include the files necessary for running the code you write to perform the tasks described here. Do not include any files from other assignments (e.g., your Unix or R assignments) or from the course tutorials. 
* Download the necessary files from the [`Python_Assignment`](https://github.com/EEOB-BioData/BCB546_Spring2025/tree/main/assignments/Python_Assignment) folder in the course repository and add them to your own repository where you will submit this assignment. The files required are:
    * `sequence_translate.py`
    * `penguins_mass.csv`
    * `penguins_cytb.fasta`
* Be sure to install the Biopython library so that this script will work. 
* You will find explicit instructions and the code/comments from Dr. X in the file called `sequence_translate.py`. Specific tasks are enumerated (1-11) in the script file with comments at the top as well as in the necessary places in the code.
* Create a new Jupyter notebook in your repository. 
    * In this notebook you will carefully annotate and execute the code already provided in `sequence_translate.py`. If you're unfamiliar with a bit of code, you can find lots of resources and information online. Be sure to cite information appropriately (by providing URLs and other relevant references). 
	* Finding new functions via web or message board searching is preferred.
    * You must also write the missing code that is outlined by pseudocode and comments.
    * Remember to document everything very clearly, following documentation guidelines we covered in class.
* Commit and push your completed Python analyses in the Jupyter notebook to your repository
    * Be certain that your GitHub repository is self contained--i.e., if someone cloned it they could run your whole Jupyter notebook without any trouble.
* Submit the URL to your git repository to Python Assignment on Canvas by the end of the day on May 7, 2025.



## Running tests

Use `pytest -q` from the repository root to run the unit tests.
