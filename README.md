# Calibration of a Local Volatility Surface

## Project Overview 
This research project explores the calibration of a local volatility surface, with a specific focus on the SSVI (Surface Stochastic Volatility Inspired) parameterisation. An emphasis is placed on the advantages of using a local volatility model that transforms the traditionally constant volatility of the Black-Scholes model into a dynamic function of time and stock price, which aims to better align with observed market prices across different maturities and strikes. This provides a more accurate representation of financial markets.

The primary objectives of this project include investigating the construction of an implied volatility surface that accurately captures the skew and term structure observed in the market, while ensuring adherence to arbitrage-free conditions. The model is calibrated to S&P500 data to estimate key parameters and is tested through the pricing of vanilla European options to assess its ability to recover market values.

This involves the tasks of:
* Preparing the market data for calibration
* Calibrating the model(/s) to the data
* Analysing the performance of the calibration
* Testing the recovery of market prices

## Code Overview
Although this code base might appear maze-like with its structure of folders, sub-folders, and various connected files, the intention was to: 1) Separate the tasks necessary to complete the project into steps which can be analysed individually, and 2) Set up the files in a way that the method can be repeated for different data sets and models by just changing a few variable names. 

The method follows the sequence:
1. Prepare the data for calibration (*Data_prep/*)
  - Run *SPX_analysis*
  - Run *SPX_bid_ask_spread*
  - The results (which will be used in the next steps) are then accessible from the *Data_prep/Data* folder
2. Calibrate the model(/s) (*Calibration/*)
  - Run the various files for each model variation
  - The results are then accessible from the *Calibration/Calibration_results* folder
3. Analyse the model performance (base folder)
  - Run the *ImpliedVolatilities* file to compute the values and errors for a particular model variation
  - Run *ModelComparisons* file to prepare visualisations that compare the different model variations' performance
4. Test the recovery of market prices (*Pricing/*)
  - Run the *EuropeanOptionPricing* file
