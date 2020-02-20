# glycans_RSSSO2015

# Ultimate IgG test

Goal: Locating sources of variation, exploring experimental variation, benchmark set for batch correction and normalization procedures development, exploring the stability of IgG 

Experiment description: Blood samples have been taken from 10 people (5 males and 5 females, approximately matched by age) during three time periods (day 1, day 2 and day 9). From these samples 3 plates were prepared by three different lab analysts (labPerson) in month 1 (labPeriod==1) and 3 more plates were prepared by the same lab analysts a month later (labPeriod==2). Therefore, accounting for variability coming from different analysts preparing the samples and for the time between the preparation of plates during an experiment with many plates. All plates had the same (randomized layout) containing IgG samples in triplicates and additional 14 standards (Genos internal standards). Hence, the analysts had to prepare the same plates.

Comments: Unfortunately, the designed turned out to be unbalanced since for three people blood samples from one day are missing. However, this should only affect the longitudinal study of the stability of IgG. Description of the data (tables):

File “ultimate_igg-areas-20150417.xlsx”
◦ Plate – plate identifier

◦ labPerson – identifier for the analysts preparing a plate

◦ labPeriod – identifier for the month (time) of plate prepration

◦ gid – Genos id for the analyzed samples

◦ Low intensity (by feeling) – denotes if the analysts doing the QC on chromatograms noticed that the intensities for a sample were smaller than he is accustomed to see – these samples pass the QC

◦ Weird chromatogram - denotes if the analysts doing the QC on chromatograms noticed that the chromatograms for a sample were different than he is accustomed to see – these samples can pass the first step of the QC but the person doing the later (statistical) analysis on these samples should be aware that the problem exists and should discarded them if they&#39;re making problems during the analysis

◦ Splitted peaks/very bad chr – denotes if the analysts doing the QC on
chromatograms finds the results for a sample really bad – these samples are not fit
to pass the QC

File “ultimate_igg-phenotypes-20150417.xlsx” - contains sex and age phenotypes for the people whose blood was sampled.

File “ultimate_igg-sample_description-20150417.xlsx” - translates gid (Genos id) to Person/Day info
