# Virus-persistence-after-spillover
Data and code for manuscript: Early epidemiological characteristics explain the chance of persistence following virus spillover events.

The R File SpilloverCharacteristics.DataPrep5.14.2025 contains code to organize and graph the raw data. It contains the code to calculate maximum likelihood prevalence and shedding from the raw qPCR and shedding results. In addition, it contains code to calculate the median corrected infection intensity.

Use data in spill.char.dataset1.csv to run the code PredictingPersistenceCorrelativeModel.R. This code runs the statistics for the correlative models described in the manuscript. 

Use data in spill.char.dataset2.csv to run the code PredictingPersistenceMechModel.R. This code runs the mechanistic model and tests its quality in a suite of models. 
