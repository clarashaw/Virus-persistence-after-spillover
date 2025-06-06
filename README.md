# Virus-persistence-after-spillover
Data and code for manuscript: Early epidemiological characteristics explain the chance of persistence following virus spillover events.

The R File SpilloverCharacteristics.DataPrep5.14.2025 contains code to clean, organize, and graph the raw data. It contains the code to calculate maximum likelihood prevalence and shedding from the raw qPCR and shedding results (L163-224). In addition, it contains code to calculate the median corrected infection intensity (L253-295). With this document, you will need to use the data files:
FCXexperiment.population.csv
FCYexperiment.population.csv
FCZexperiment.population.csv
FDIexperiment.population.csv
FCX.strips.csv
FCY.strips.csv
FCZ.strips.csv
FDI.strips.csv
FCX.shedding.csv
FCY.shedding.csv
FCZ.shedding.csv
FDI.shedding.csv
Strainsandblocks.csv
TCID50.csv

The R file PredictingPersistenceCorrelativeModel.R contains the code for the correlative models described in the manuscript. SpilloverCharacteristics.DataPrep5.14.2025 will make the .csv file to use with this code, but I have also uploaded it separately (spill.char.dataset1.csv)

The R file PredictingPersistenceMechModel.R contains the code for the mechanistic model and the suite of models used to test its quality and understand where the model falls short. SpilloverCharacteristics.DataPrep5.14.2025 will make the .csv file to use with this code, but I have also uploaded it separately (spill.char.dataset2.csv)
