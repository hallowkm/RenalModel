library(RxODE)

source("modelfile.R")
load("model_struct.saved")

m1 <- RxODE(model = ode, modName = "mod1")

#load basecase parameters
source("calcNomParams_human.R")
theta=calcNomParams()

###Initial conditions - do NOT change order!!!
#Order must match order in model file
#labels are not used by RxODE to match init to compartment
inits <- c(AngI=8.164, AngII=5.17,
AT1_bound_AngII=16.6, AT2_bound_AngII = 5.5, plasma_renin_concentration=17.845,
blood_volume_L = theta["blood_volume_nom"],
 extracellular_fluid_volume=theta["ECF_nom"],
sodium_amount= as.numeric(theta["blood_volume_nom"])*as.numeric(theta["ref_Na_concentration"]), 
ECF_sodium_amount= as.numeric(theta["ECF_nom"])*as.numeric(theta["ref_Na_concentration"]), 
tubulo_glomerular_feedback_effect=1,
normalized_aldosterone_level_delayed=1, 
preafferent_pressure_autoreg_signal=1, 
glomerular_pressure_autoreg_signal=1, 
cardiac_output_delayed=theta["CO_nom"],
CO_error=0, Na_concentration_error = 0, 
normalized_vasopressin_concentration_delayed = 1,
F0_TGF=theta["nom_LoH_Na_outflow"], 
P_bowmans=theta["Pc_pt_mmHg"], 
oncotic_pressure_difference=theta["nom_oncotic_pressure_difference"],
SN_macula_densa_Na_flow_delayed = as.numeric(theta["nom_LoH_Na_outflow"])/as.numeric(theta["baseline_nephrons"]),
serum_creatinine = as.numeric(theta["equilibrium_serum_creatinine"])*as.numeric(theta["blood_volume_nom"]), cumNaExcretion = 0, cumWaterExcretion = 0, cumCreatinineExcretion = 0
)

dev.new()

#calculate steady state conditions
ev1 <- eventTable()
ev1$add.sampling(seq(0,100000,by=10))
x <- m1$run(theta, ev1, inits=inits)

#Store final values as new initial conditions
inits=x[dim(x)[1],2:(length(inits)+1)]

#Run model at steady state
x <- m1$run(theta, ev1, inits=inits)

#Check that model is at steady-state
dev.new()
par(mfrow=c(1,2))
matplot(x[,"mean_arterial_pressure_MAP"],type="l")
plot(x[,"GFR_ml_min"],type="l")



