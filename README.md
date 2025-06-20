This zip file contains all R code required to reproduce the figures and tables of this article. While code to analyze the data used in the data analysis section is provided, the data itself is not provided. Below are descriptions of each of the files in this folder.

## options.csv
List of options for simulation parameters.

## X_pregen.R
Generate data for simulations and save, along with a chosen lambda_max which will be used to generate ROC curves.

## main.R
Main script which will run simulations on pre-generated data.

## meMethods.R
Helper functions for generating error-prone data, and for correcting for measurement error bias in the variance-covariance matrix.

## cdlmmLasso.R
Set of functions to implement the coordinate-descent LMM Lasso algorithm from the paper, as well as the algorithm which selects the BIC-optimal model automatically.

## helperfunctions.R
Contains a function which returns a covariance matrix with specified structure.

## simfunctions.R
Functions to run simulations, either in parallel or serially.

## make_res.R
Makes a table and a figure to show estimation error and ROC performance from simulation study.

## data_parsing.R
Code to preprocess IDATA study data

## lmmLasso_analysis.R
Code to analyze preprocessed IDATA study data and produce the results in the paper.