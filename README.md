# SWG
Title: Event Attribution of Climate Changes with Dynamically driven-Stochastic Weather Generators

## How to run the code
Please download and RUN `slp_SD_diagnosis.R`
- [slp_SD_diagnosis.R](https://github.com/cfcforever/SWG/blob/master/slp_SD_diagnosis.R)
- dossier [**function**](https://github.com/cfcforever/SWG/tree/master/function)
- dossier [**DATA**](https://github.com/cfcforever/SWG/tree/master/DATA)

## How to construct conditional SWG
We use **tmean** as example (see [slp_SD_diagnosis](https://github.com/cfcforever/SWG/blob/master/slp_SD_diagnosis.R)), 
1. to select the predictors (e.g., using the dynamical system on sea level pressure over the Northern Atlantic)
2. to estimate the parameters for mean and sd conditional to the dynamical system
3. to simulate the mean and sd conditional to the dynamical system
4. to generate the simulations by using the simulated mean and sd

5. to repeat step 2-4 to generate the simualtions of non-conditional model
> We comapre the conditional model to the non-conditional model

## How to use conditaional SWG for event attribution of climate changes
1. FAR
2. Changes of intensity
3. Anomaly maps


## To complete
### Coordinations

| Name                      | Domain                |
| ------------------------- |:---------------------:|
| Northern Atlantic         | [22N, 71N]x[81W, 51E] |
| Northern Atlantic (small) | [30N, 71N]x[30W, 51E] |
| Europe                    | [35N, 60N]x[10W, 25E] |


### Function
- [fun_estimation](https://github.com/cfcforever/SWG/blob/master/function/fun_estimation_t2m.R)
- [fun_simulation](https://github.com/cfcforever/SWG/blob/master/function/fun_simulation_t2m.R)
- [plot_worldmap](https://github.com/cfcforever/SWG/blob/master/function/plot_worldmap.R)

### Anomaly 
