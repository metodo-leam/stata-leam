# metodo-leam/stata - Stata user-written commands by the Laboratori d'Estadística Aplicada i Modelització (UAB)

## Content

The repository `stata` provides Stata user-written commands by the Laboratori d'Estadística Aplicada i Modelització, Universitat Autònoma de Barcelona (UAB). 

| Command | Description |
| --- | --- |
| `agree` | Agreement: Bland-Altman & Passing-Bablok methods |
| `allsets` | All Possible Subsets: Linear, Logistic, Cox proportional hazards and Competing-risk (Fine-Gray) regression |
| `chisqi` | Goodness of fit Chi-squared test |
| `cir` | Confidence Intervals for Pearson and Spearman Correlation | 
| `cohenkap` | Kappa and Weighted kappa |
| `confound` | Modelling confounding in Linear, Logistic, Cox proportional hazards and Competing-risk (Fine-Gray) regression |
| `dc` | Data Check |
| `dcreport` | Data Check - Incidence Report |
| `dt` | Diagnostic Tests |
| `dtroc` | ROC Analysis & Optimal Cutoff Point |
| `intcon` | Internal Consistency |
| `mar` | Meta-Analysis: OR,RR,RD,IR,ID,B,MD,R combined |
| `nsize` | Sample Size & Power |
| `pwkwallis` | Kruskal-Wallis equality-of-populations rank test & pairwise comparisons |
| `rndseq` | Generation of Random Sequences |
| `rtrend` | Trend test (2xk or kx2 table) for frequency and person-time data |
| `scct` | Stochastic Curtailment: Clinical trials |
| `sta` | Association measures (frequency, person-time & paired data) |
| `statmis` | Statistics of missing values |

Each command includes a help file and a dialog box.

These commands have been tested with Stata versions 12 to 19.

### Install

To install one command execute the syntax below, click on the name an follow the instructions.

```stata
net from "https://raw.githubusercontent.com/metodo-leam/stata/master"
```

You may install using `net install`. For example, to install the `agree` command execute:

```stata
net install agree, from("https://raw.githubusercontent.com/metodo-leam/stata/master") replace
```

To install other commands, replace `agree` with the name of the desired command.

### Background

These commands were developed as part of the studies _Metodología de la investagación: Diseño y Estadística en Ciencias de la Salud_, at the Laboratori d'Estadística Aplicada i Modelització (UAB).

## Citation

Please use the Vancouver reference suggested in the _Authors_ section of the `help` of each command.

### Contact

[metodo.campus@gmail.com](mailto:metodo.campus@gmail.com)
