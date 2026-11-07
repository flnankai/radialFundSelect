# Replicating Figure 4 and the U.S. equity fund results

The application uses the panel distributed with Jiti Gao, Bin Peng and Yayi Yan (2024), “Higher-Order Expansions and Inference for Panel Data Models,” *Journal of the American Statistical Association* 119, 2760–2771. They constructed the panel from the CRSP Survivor-Bias-Free Mutual Fund Database. The authors provide the research archive on their publication page. We found no license in the archive granting permission to redistribute its third-party fund observations. Obtain your own authorized copy and keep it outside the package and repository. The package contains no individual fund returns or downloaded French data. [Paper](https://doi.org/10.1080/01621459.2023.2277411) · [Author's publication page, “Matlab Code” next to this title](https://sites.google.com/view/yayiyan/yayi-yan) · [Author's linked archive](https://drive.google.com/file/d/1-Frlu2Ac1RYfCAuh2ZK_Nivm-tSDMU7u/view)

1. Unpack the authors' archive locally and locate `CaseStudy/mutual_fund.mat`. The archive directory may be nested inside a release folder. The package reads this MAT file directly in R. Its matrices `R` and `F` contain 668 months and have 726 and four columns respectively. Both are already in decimal units. The replication uses rows 301 through 668, for February 1987 through September 2017. The four columns of `F` are market excess return, SMB, HML and momentum.
2. Download the **July 2024 CSV archive** under historical U.S. Fama/French 3 Factors. The paper used the RF column from the June 2024 data cut. The July 2024 archive's 368 monthly RF observations over the study period match the paper's inputs exactly. You can pass the downloaded ZIP directly. The reader extracts only the RF column, aligns it by YYYYMM, and converts percentages to decimals. It retains the four factors from the authors' MAT file. [French historical archives](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/f-f_factors_archive.html) · [Official July 2024 ZIP](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/Historical_Archives/08%202024%20Update/ftp/F-F_Research_Data_Factors_CSV.zip)
3. Run the following commands using your local file paths. The output directory must already exist. Install the optional MAT reader once with `install.packages("R.matlab")`.

```r
library(radialFundSelect)
emp <- reproduce_empirical(
  data = "path/to/CaseStudy/mutual_fund.mat",
  rf_file = "path/to/F-F_Research_Data_Factors_CSV.zip",
  figure = "empirical_T69.pdf"
)
emp$summary
```

If you already have the paper's three aligned, headerless CSV inputs, place them together in a private directory and run `reproduce_empirical("path/to/private/directory", figure = "empirical_T69.pdf")`. This route uses base R and the installed package. The files are `fund.csv` (726 funds by 368 calendar months, decimal gross-of-risk-free returns reported net of fund fees and expenses), `factor.csv` (368 months by four decimal factors in the order above), and `rf.csv` (368 monthly risk-free percentages, one per line). Preserve the source's zero codes in the fund panel. The documentation for `reproduce_empirical` describes how these are handled separately during training and investment.

The result contains one row per holding frequency and method in `summary`, formation counts in `periods`, and aggregate monthly paths in `wealth` and `stress_wealth`. Rates are decimals and drawdowns are positive losses. The 3-, 6- and 9-month designs have 99, 49 and 33 complete holding periods. IN-BH's compound annual growth rates round to 4.96%, 4.00% and 3.47%. Its maximum drawdowns are 5.44%, 5.35% and 13.84%. The common-endpoint columns value all holdings through April 2017. Run `plot_empirical_wealth(emp, "empirical_T69.pdf")` to redraw the figure using the computed paths.
