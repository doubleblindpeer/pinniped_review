## R files descriptions
pinniped_map_temporal_figs.R - Plots the world map and temporal observation coverage of the real pinniped database.

simulation_distance_covariance.R - Generates the simulated data for the simulation experiment. Also calculates spatial distances among intraspecific populations.

MARSS_model_estimation_sim_real.R - Estimates species trends and state-space parameters for simulated data and real dataset. Outputs the species trend plot(s).

sim_processing.R - Uses model outputs of simulated data to evaluate the model selection framework and parameter estimation accuracy.

lambda_calculations.R - Calculates and plots the total species trend lambdas (annual rates of change). Associated with intra- and inter-specific lambda plot, trend in annual growth rates of all pinniped species plot, and IUCN plot.

compare_model_types_species_trend.R - Compares species trends of the real database to their LPI counterparts.

## Data files descriptions
Pinniped database for publication.csv - Real database of 555 time series for pinniped abundances.

seal_IUCN_gen_length.xlsx & pinniped_generation_lengths.xlsx - IUCN characteristics and pinniped generation lengths data (Source: IUCN, https://natureconservation.pensoft.net/article/1343/)

## Setting up the MARSS model
Download CmdStanR and follow the help guide: https://mc-stan.org/cmdstanr/articles/cmdstanr.html

Stan files are located in the cmd_stan_files folder
