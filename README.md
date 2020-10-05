# worm_survival
This R script produces survival stats and a Kaplan-Meier plot (originally created for wax worm experiments). 
<p>&nbsp;</p>

<b>Usage</b>
```
Rscript worm_survival.R WORM_COUNT_DATA.csv
```
<p>&nbsp;</p>

<b>Input</b>
- Survival Data CSV file - Column 1 is the day starting with day 0 (the number of worms at start of experiment). Any additional columns wiill be counted as individual test conditions. You can input as many test conditions as you want. An example ("worms1.csv") is located in the test folder. 
<p>&nbsp;</p>

<b>Output</b>
filename | description
-------- | -------------
LT50.csv | CSV contaning the LT50 (time, in days, when 50% of worms are dead) for all test conditions
LogRankPvalues.csv | CSV contaning the log rank p-values between each pair of test conditions (a test for significance) 
cox_plot.jpeg | a plot containing hazard ratios, using the first test condition as the negative control reference
survival_plot.jpeg | Kaplan-Meier plot of survivial probability over time
FormattedSurvivalData.csv | This is the survival data reformatted so the Kaplan-Meier plot can be produced; fustat=0 if the worm survived the whole study, fustat=1 if the worm died before the end of the study; futime=day when worm died or day when the study ended (if the worm survived the whole time)
Rplots.pdf | empty file that needs to be deleted
<p>&nbsp;</p>

* Kaplan-Meier Survival Plot

<img src="https://github.com/amcrabtree/worm_survival/blob/main/images/survival_plot.jpeg" alt="drawing" width="500"/>

* Cox Hazard Plot

<img src="https://github.com/amcrabtree/worm_survival/blob/main/images/cox_plot.jpeg" alt="drawing" width="500"/>
