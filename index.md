
# DRpower

The *DRpower* package is designed to assist with the development of *Plasmodium falciparum* *pfhrp2/3* gene deletion studies. It can be used in the following ways:

- **In the design phase** to choose a number of sites and sample size per site to achieve a desired level of power. This is a critical step that ensures the conclusions of a study are statistically valid
- **In the analysis phase** to perform calculations to estimate the prevalence of *pfhrp2/3* deletions and to test whether the prevalence is over a set threshold (5% by default) at the domain level

<p style="text-align:center;">
  <img src="https://raw.githubusercontent.com/mrc-ide/DRpower/master/R_ignore/outputs/homepage_figure.png" height="338px" width="525px" />
  <br>
  Example output from *DRpower*. This distribution gives the probability of pfhrp2/3 gene deletions being at any given prevalence in the domain. The red shaded region represents prevalence being above the 5% threshold - in this case there is an 87.4% chance of being above the threshold.
</p>

The simplest way to interact with the package is through the [interactive web app](https://shiny.dide.ic.ac.uk/DRpower-app/), which contains all the main functions of the package in a user-friendly environment. If you are looking for a simple way to design a study or analyse your data then you should start with the app. Working directly with this package is recommended for people who want to do more advanced analyses, for example considering multiple end-points, or who want to set up reproducible pipelines.

The ["How it works"](articles/rationale1_background.html) section gives a thorough explanation of the rationale behind the *DRpower* method. If you want to gain a deeper understanding of the statistical issues that motivate this software then you should start here.

For those wanting to jump straight into using the software, see the [installation instructions](articles/installation.html) followed by a series of [tutorials](articles/tutorial_design.html) on how to perform various analyses.
