# worm_survival
This R script produces a PDF report of survival stats and plots, as well as a CSV containing the data reformatted for survival analysis. This script was originally created for evaluating toxicity of substances administered to wax moth larvae. 
<p>&nbsp;</p>

<b>Usage (unix/bash interface)</b>
>Make sure that the rmarkdown library is installed on your R CLI. If in doubt, run the following command: `Rscript -e 'install.packages("rmarkdown")'`. The worm_survival.Rmd script will automatically download any other necessary packages for its execution. 
```
Rscript -e 'rmarkdown::render("worm_survival.Rmd", params=list(data="test/worms1.csv"))'
```
<p>&nbsp;</p>

<b>Input</b>
- Survival Data CSV file - Column 1 is the day starting with day 0 (the number of worms at start of experiment). Any additional columns wiill be counted as individual test conditions. You can input as many test conditions as you want. An example ("worms1.csv") is located in the test folder. 
<p>&nbsp;</p>

<b>Output</b>
filename | description
-------- | -------------
FormattedSurvivalData.csv | This is the survival data reformatted so the Kaplan-Meier plot can be produced; fustat=0 if the worm survived the whole study, fustat=1 if the worm died before the end of the study; futime=day when worm died or day when the study ended (if the worm survived the whole time)
worm_survival.pdf | This is a PDF report containing survival plots and other statistics. 
<p>&nbsp;</p>

* Example of a Kaplan-Meier Survival Plot

<img src="https://github.com/amcrabtree/worm_survival/blob/main/images/survival_plot.jpeg" alt="drawing" width="500"/>

* Example of a Cox Hazard Plot

<img src="https://github.com/amcrabtree/worm_survival/blob/main/images/cox_plot.jpeg" alt="drawing" width="500"/>
