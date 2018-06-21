calcNomParams <- function(){


########################################################################################
#Parameters of normal human physiology based on literature and commmon medical knowledge
########################################################################################

####Systemic parameters
nominal_map_setpoint=93  		#mmHg
CO_nom= 5					#L/min
ECF_nom = 15				#L
blood_volume_nom = 5			#L
Na_intake_rate=100/24/60		#mEq/min  - 100mmol/day or 2300 mg/day
nom_water_intake = 2.1			#L/day
ref_Na_concentration=140 	   	#mEq/L
plasma_protein_concentration = 7   	#g/dl
equilibrium_serum_creatinine=0.92	#mg/dl
P_venous=4					#mmHg
R_venous=3.4				#mmHg
nom_right_atrial_pressure=0.87 	#mmHg
nom_mean_filling_pressure=7 		#mmHg
venous_compliance = 0.13


####Renal parameters
nom_renal_blood_flow_L_min = 1	#L/min
baseline_nephrons=2e6
nom_Kf=3.9					#nl/min*mmHg
nom_oncotic_pressure_difference= 28 #mmHg
P_renal_vein=4  				#mmHg

#Renal Vasculature
nom_preafferent_arteriole_resistance= 19 	#mmHg
nom_afferent_diameter=1.5e-5		 	#mmHg
nom_efferent_diameter=1.1e-05 		#mmHg

#Renal Tubules
Dc_pt_nom  = 27e-6			#m
Dc_lh = 17e-6				#m
Dc_dt = 17e-6				#m
Dc_cd = 22e-6				#m

L_pt_s1_nom = 0.005			#m
L_pt_s2_nom = 0.005			#m
L_pt_s3_nom =0.004			#m
L_lh_des = 0.01 				#m
L_lh_asc = 0.01 				#m
L_dct = 0.005				#m	
L_cd = L_lh_des	

tubular_compliance = 0.2	
Pc_pt_mmHg = 14				#mmHg
Pc_lh_des_mmHg = 10.5			#mmHg
Pc_lh_asc_mmHg = 7			#mmHg
Pc_dt_mmHg = 3				#mmHg
Pc_cd_mmHg = 2				#mmHg
P_interstitial_mmHg = 5
nominal_pt_na_reabsorption=0.7	#fraction
nominal_loh_na_reabsorption = 0.8	#fraction
nominal_dt_na_reabsorption=0.5	#fraction
LoH_flow_dependence = 1


####RAAS Pathway parameters
concentration_to_renin_activity_conversion_plasma = 61 
nominal_equilibrium_PRA = 1000 	 	#fmol/ml/hr
nominal_equilibrium_AngI = 7.5 		#fmol/ml
nominal_equilibrium_AngII = 4.75 		#fmol/ml
nominal_renin_half_life = 0.1733		# (hr)
nominal_AngI_half_life = 0.5/60 		#(hr)
nominal_AngII_half_life = 0.66/60 		#(hr)
nominal_AT1_bound_AngII_half_life = 12/60 #hr
nominal_AT2_bound_AngII_half_life = 12/60 #hr
ACE_chymase_fraction = 0.95     		#% of AngI converted by ACE. The rest is converted by chymase
fraction_AT1_bound_AngII = 0.75    		#assume AngII preferentially binds to AT1 vs AT2


########################################################################################
#The following parameters are calculated at equilibrium using the parameters above
########################################################################################

#This pressure is the setpoint that determines the myogenic response of the preafferent vasculature
nom_preafferent_pressure = nominal_map_setpoint - nom_renal_blood_flow_L_min*nom_preafferent_arteriole_resistance;

#This pressure is the setpoint that determines the myogenic response of the afferent vasculature
nom_glomerular_pressure = nom_preafferent_pressure - nom_renal_blood_flow_L_min*(L_m3*viscosity_length_constant/(nom_afferent_diameter^4)/baseline_nephrons);

#This pressure is the setpoint that determines the tubular pressure-natriuresis response 
nom_postglomerular_pressure = nom_preafferent_pressure - nom_renal_blood_flow_L_min*(L_m3*viscosity_length_constant*(1/(nom_afferent_diameter^4)+1/(nom_efferent_diameter^4))/baseline_nephrons);

RIHP0 = nom_postglomerular_pressure 	

# The rate of sodium excretion must equal the rate of sodium intake. Sodium reabsorption rates vary along the tubule, but based on literature
# measurements we have a good, and literature data provides estimates for these rates. However, there is a precise
# rate of sodium reabsorption required to achieve the equilibrium defined by the parameters above.
# Assuming that reabsorption rates are known in all but one segment of the tubule, the exact rate
# of reabsorption of the remaining segment can be calculated. We chose to calculate the CD rate of reabsorpion based on estimates for
# PT, LoH, and DT reabsorption.  
nom_GFR = nom_Kf*(nom_glomerular_pressure - nom_oncotic_pressure_difference - (Pc_pt_mmHg+P_interstitial_mmHg))/nL_mL*baseline_nephrons;
nom_filtered_sodium_load = nom_GFR/L_mL*ref_Na_concentration;
nom_PT_Na_outflow = nom_filtered_sodium_load*(1-nominal_pt_na_reabsorption);

nom_Na_in_AscLoH = nom_PT_Na_outflow/baseline_nephrons;
AscLoH_Reab_Rate =(2*nominal_loh_na_reabsorption*nom_Na_in_AscLoH)/L_lh_des; #osmoles reabsorbed per unit length per minute. factor of 2 because osmoles = 2

nom_LoH_Na_outflow = nom_PT_Na_outflow*(1-nominal_loh_na_reabsorption);
nom_DT_Na_outflow = nom_LoH_Na_outflow*(1-nominal_dt_na_reabsorption);
nominal_cd_na_reabsorption = 1-Na_intake_rate/nom_DT_Na_outflow;



#RBF = (MAP - P_venous)/RVR. Given MAP, P_venous, RBF, and preafferent, afferent, and efferent resistances, the remaining peritubular resistance at steady state can be determined
nom_RVR = (nominal_map_setpoint - P_venous)/nom_renal_blood_flow_L_min
nom_peritubular_resistance = nom_RVR - (nom_preafferent_arteriole_resistance + L_m3*viscosity_length_constant*(1/nom_afferent_diameter^4+1/nom_efferent_diameter^4)/baseline_nephrons);


#Given the values for baseline MAP and CO above, the baseline TPR required to maintain this MAP and CO can be calculated. Since TPR includes renal vascular resistance, the baseline systemic (non-renal) resistance
#can be calculated from this TPR and the values for baseline renal resistances defined above. 
nom_TPR = nominal_map_setpoint/CO_nom
nom_systemic_arterial_resistance= nom_TPR-R_venous

#Creatinine synthesisrate at equilibrium
creatinine_synthesis_rate  = equilibrium_serum_creatinine * dl_ml * nom_GFR #Units: mg/min


####RAAS Pathway parameters
#Values for half lives and equilibrium concentrations of RAAS peptides available in the literature and 
# defined above to calculate nominal values for other RAAS parameters not available in the literature:
#ACE activity
#Chymase activity
#AT1 receptor binding rate
#AT2 receptor binding rate
#equilibrium AT1_bound_AngII
#These values are then assumed to be fixed unless specified otherwise.
#Calculating these nominal parameter values initially in a separate file is required so that these parameters can then be varied independently in the main model
nominal_equilibrium_PRC = nominal_equilibrium_PRA/concentration_to_renin_activity_conversion_plasma
nominal_AngI_degradation_rate = log(2)/nominal_AngI_half_life #/hr
nominal_AngII_degradation_rate = log(2)/nominal_AngII_half_life #/hr
nominal_AT1_bound_AngII_degradation_rate = log(2)/nominal_AT1_bound_AngII_half_life
nominal_AT2_bound_AngII_degradation_rate = log(2)/nominal_AT2_bound_AngII_half_life
#ACE converts 95% of AngI, chymase converts the rest
nominal_ACE_activity = (ACE_chymase_fraction*(nominal_equilibrium_PRA - nominal_AngI_degradation_rate*nominal_equilibrium_AngI)/nominal_equilibrium_AngI)#Therapy_effect_on_ACE
nominal_chymase_activity = (1-ACE_chymase_fraction)*(nominal_equilibrium_PRA - nominal_AngI_degradation_rate*nominal_equilibrium_AngI)/nominal_equilibrium_AngI
#75% of bound AngII is AT1, the rest is AT2
nominal_AT1_receptor_binding_rate = fraction_AT1_bound_AngII*(nominal_equilibrium_AngI*(nominal_ACE_activity+nominal_chymase_activity)-nominal_AngII_degradation_rate*nominal_equilibrium_AngII)/nominal_equilibrium_AngII
nominal_AT2_receptor_binding_rate = (1-fraction_AT1_bound_AngII)*(nominal_equilibrium_AngI*(nominal_ACE_activity+nominal_chymase_activity)-nominal_AngII_degradation_rate*nominal_equilibrium_AngII)/nominal_equilibrium_AngII
nominal_equilibrium_AT1_bound_AngII = nominal_equilibrium_AngII*nominal_AT1_receptor_binding_rate/nominal_AT1_bound_AngII_degradation_rate
nominal_equilibrium_AT2_bound_AngII = nominal_equilibrium_AngII*nominal_AT2_receptor_binding_rate/nominal_AT2_bound_AngII_degradation_rate

########################################################################################
#The following parameters were determined indirectly from many different literature studies on the response
#various changes in the system (e.g. drug treatments, infusions of peptide, fluid, sodium, etc.....)
########################################################################################

#Effects of AT1-bound AngII on preafferent, afferent, and efferent resistance, and aldosterone secretion
AT1_svr_slope = 0
AT1_preaff_scale = 0.5
AT1_preaff_slope = 7 
AT1_aff_scale=0.5
AT1_aff_slope=7
AT1_eff_scale=0.3
AT1_eff_slope=7
AT1_PT_scale = 0.1
AT1_PT_slope = 7
AT1_aldo_slope = 0.05


#Effects of Aldosterone on distal and collecting duct sodium reabsorption
nominal_aldosterone_concentration=85
aldo_DCT_scale=0
aldo_DCT_slope = 0.5
aldo_CD_scale=0.3
aldo_CD_slope = 0.5
aldo_renin_slope =-0.05

#Na and water transfer between blood, ECF
Q_water = 1
Q_Na = 1

#Osmolarity control of vasopressin secretion
Na_controller_gain=2  
Kp_VP = 0.05
Ki_VP = 0.00002

nom_ADH_urea_permeability = .98
nom_ADH_water_permeability = .98

#Effects of Vasopressin on water intake and reabsorption
nominal_vasopressin_conc=4
water_intake_vasopressin_scale = 0#1.5
water_intake_vasopressin_slope = -0.5


#Magnitude and Steepness of tubuloglomerular feedback
S_tubulo_glomerular_feedback=0.7
F_md_scale_tubulo_glomerular_feedback=6
MD_Na_concentration_setpoint = 62.4

#Effect of macula densa sodium flow on renin secretion 
md_renin_A = 1
md_renin_tau = 2

#Responsiveness of renal vasculature to regulatory signals
preaff_diameter_range=0.25
afferent_diameter_range=1.2e-05 
efferent_diameter_range=3e-06 
preaff_signal_nonlin_scale=3
afferent_signal_nonlin_scale=3
efferent_signal_nonlin_scale=3

#RAAS pathway (these parameters can be set to different values than used to calculate the equilibrium state above)
AngI_half_life=0.008333 
AngII_half_life=0.011 
AT1_bound_AngII_half_life=0.2 
AT1_PRC_slope=-1.2 
AT1_PRC_yint=0
AT2_bound_AngII_half_life=0.2
concentration_to_renin_activity_conversion_plasma=61
fraction_AT1_bound_AngII=0.75
nominal_ACE_activity=48.9
 nominal_AT1_receptor_binding_rate=12.1
nominal_AT2_receptor_binding_rate=4.0 
nominal_chymase_activity=1.25  
nominal_equilibrium_AT1_bound_AngII=16.63
nominal_equilibrium_PRC=16.4 
renin_half_life=0.1733 

#Transfer constants for ODEs - determine speed of processes
C_aldo_secretion=1000
C_P_bowmans = 1000
C_P_oncotic = 1000

C_tgf_reset=0
C_cardiac_output_delayed=.001
C_co_error=0.00001
C_vasopressin_delay = 1

C_md_flow = 0.001 #Time delay between MD sodium flow and renin secretion
C_tgf=1
C_na_excretion_na_amount=-1
C_na_intake_na_amount=1
C_urine_flow_ecf_volume=-1
C_water_intake_ecf_volume=1
C_Na_error=1/6
C_serum_creatinine = 1



#Therapy effects
HCTZ_effect_on_DT_Na_reabs = 1 
HCTZ_effect_on_renin_secretion = 1
DRI_effect_on_PRA = 1
CCB_effect_on_preafferent_resistance = 1
CCB_effect_on_afferent_resistance = 1
CCB_effect_on_efferent_resistance = 1
MR_antagonist_effect_on_aldo_MR = 1
pct_target_inhibition_ARB = 0 
pct_target_inhibition_ACEi = 0 


#Metabololic tissue autoregulation of cardiac output
tissue_autoreg_scale=1
Kp_CO=1.5
Ki_CO=30

 
#Normalized_aldo_secretion
K_Na_ratio_effect_on_aldo = 1; 

#Renal autoregulation of glomerular pressure
gp_autoreg_scale=0
preaff_autoreg_scale = 0.5
myogenic_steepness=2


#Pressure natiuresis effect 

max_pt_reabs_rate = 0.995
pressure_natriuresis_PT_scale = 3
pressure_natriuresis_PT_slope = 1

pressure_natriuresis_LoH_scale = 3
pressure_natriuresis_LoH_slope = 1

pressure_natriuresis_DCT_scale = 3
pressure_natriuresis_DCT_slope = 1

max_cd_reabs_rate = 0.995
pressure_natriuresis_CD_scale = 3
pressure_natriuresis_CD_slope=1



#Constants and Unit conversions
nL_mL=1e+06
dl_ml=0.01
L_dL=10
L_mL=1000
L_m3=0.001
g_mg=0.001
ng_mg=1e-06
secs_mins=60
min_hr=60
hr_day=24
min_day=1440
MW_creatinine=113.12
Pi=3.1416
pi=3.14
viscosity_length_constant=1.5e-09 
gamma =  1.16667e-5;  #  viscosity of tubular fluid
mmHg_Nperm2_conv = 133.32


#Scaling parameters - can be used to parameterize model for other species
ECF_scale_species = 1
BV_scale_species=1
water_intake_species_scale = 1
CO_scale_species = 1


t=sort(ls())
param=sapply(t,names)
for (i in 1:length(t)){
  param[i]=get(t[i])
}
param$param=NULL

return(param)
}




