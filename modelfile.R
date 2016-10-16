##############
#Author: K. Melissa Hallow
#hallowkm@uga.edu
#University of Georgia
#Athens, GA
#
#Date: October 13, 2016
#
#Description: The model describes the key physiological processes involved in renal function 
#and its role in maintaining sodium and water homeostasis,  at the systems level based on the 
#governing  physiological and feedback mechanisms. The model is not meant to be an exhaustive 
#molecular to organ level model, but rather, to provide a backbone for further investigation. 


ode <- " 

#######################################Systemic Hemodynamics #################################################

#################  Systemic Vascular Resistance
#Systemic vascular resistance is a nominal value modulated by AngII and by a regulatory blood flow autoregulation signal

###Whole body autoregulation mechanism wherein TPR adjusts to maintain constant organ blood flow (and thus constant cardiac output)
#Modeled as Proportional-Integral controller of TPR, where the input signal is the cardiac output error signal
tissue_autoregulation_signal = max(0.1,1+tissue_autoreg_scale*((Kp_CO/CO_scale_species)*(cardiac_output_delayed - CO_nom)+(Ki_CO/CO_scale_species)*CO_error));

###Effect of the RAAS (AT1-bound AngII) on systemic vascular resistance. 
AT1_svr_int = 1 - AT1_svr_slope*nominal_equilibrium_AT1_bound_AngII;
AT1_bound_AngII_effect_on_SVR = AT1_svr_int + AT1_svr_slope * AT1_bound_AngII;

systemic_arterial_resistance = nom_systemic_arterial_resistance*tissue_autoregulation_signal*AT1_bound_AngII_effect_on_SVR;  

################# Cardiac Output and Mean Arterial Pressure
#Cardiac output is a function of blood volume and resistance to venous return

resistance_to_venous_return = ((8 * R_venous + systemic_arterial_resistance) / 31); 

mean_filling_pressure = nom_mean_filling_pressure + (blood_volume_L/BV_scale_species-blood_volume_nom)/venous_compliance;

cardiac_output =  mean_filling_pressure / resistance_to_venous_return; 

################# Cardiac Output and Mean Arterial Pressure

total_peripheral_resistance = systemic_arterial_resistance + R_venous;

mean_arterial_pressure_MAP = cardiac_output * total_peripheral_resistance;


####################################### Renal Vasculature #################################################

###AT1-bound AngII constricts the preafferent, afferent, and efferent arterioles
AT1_preaff_int = 1 - AT1_preaff_scale/2;
AT1_effect_on_preaff = AT1_preaff_int + AT1_preaff_scale/(1+exp(-(AT1_bound_AngII - nominal_equilibrium_AT1_bound_AngII)/AT1_preaff_slope));

AT1_aff_int = 1 - AT1_aff_scale/2;
AT1_effect_on_aff = AT1_aff_int + AT1_aff_scale/(1+exp(-(AT1_bound_AngII - nominal_equilibrium_AT1_bound_AngII)/AT1_aff_slope));

AT1_eff_int = 1 - AT1_eff_scale/2;
AT1_effect_on_eff = AT1_eff_int + AT1_eff_scale/(1+exp(-(AT1_bound_AngII - nominal_equilibrium_AT1_bound_AngII)/AT1_eff_slope));

#################  Preafferent Resistance
#The resistance of the arcuate, interlobular arterioles, and other vasculature prior the afferent arterioles is represented by a single resistance - the preafferent arteriole resistance
#The preafferent arterioles respond myogenically to changes in pressure, and also responds to AT1-bound AngII
#The dilation/constriction of the arterioles is limited, and thus the total combined effect of all regulators must saturate
preaff_arteriole_signal_multiplier = AT1_effect_on_preaff*(preafferent_pressure_autoreg_signal)*CCB_effect_on_preafferent_resistance;
preaff_arteriole_adjusted_signal_multiplier = (1/(1+exp(preaff_signal_nonlin_scale*(1-preaff_arteriole_signal_multiplier)))+0.5);
preafferent_arteriole_resistance = nom_preafferent_arteriole_resistance*preaff_arteriole_adjusted_signal_multiplier;

#################  Afferent Arteriole Resistance
#The afferent arteriole responses the tubuloglomerular feedback (calculated later), as well as to AT1-bound AngII. 
#It may respond myogenically as well. Some studies suggest the upstream portion responds myogenically while the distal portion responds to TGF. Thus, one could consider the 
#myogenically responsive portion as part of the preafferent resistance. 
#The dilation/constriction of the arterioles is limited, and thus the total combined effect of all regulators must saturate
nom_afferent_arteriole_resistance = L_m3*viscosity_length_constant/(nom_afferent_diameter^4);
afferent_arteriole_signal_multiplier = tubulo_glomerular_feedback_effect * AT1_effect_on_aff *glomerular_pressure_autoreg_signal*CCB_effect_on_afferent_resistance;
afferent_arteriole_adjusted_signal_multiplier = (1/(1+exp(afferent_signal_nonlin_scale*(1-afferent_arteriole_signal_multiplier)))+0.5);
afferent_arteriole_resistance = nom_afferent_arteriole_resistance*afferent_arteriole_adjusted_signal_multiplier;

#################  Efferent Arteriole Resistance
#The efferent arteriole responses to AT1-bound AngII.
#The dilation/constriction of the arterioles is limited, and thus the total combined effect of all regulators must saturate
nom_efferent_arteriole_resistance = L_m3*viscosity_length_constant/(nom_efferent_diameter^4);
efferent_arteriole_signal_multiplier = AT1_effect_on_eff *CCB_effect_on_efferent_resistance;
efferent_arteriole_adjusted_signal_multiplier = 1/(1+exp(efferent_signal_nonlin_scale*(1-efferent_arteriole_signal_multiplier)))+0.5;
efferent_arteriole_resistance = nom_efferent_arteriole_resistance*efferent_arteriole_adjusted_signal_multiplier;


#################  Peritubular Resistance
peritubular_resistance = nom_peritubular_resistance*baseline_nephrons;

#################  Renal Vascular Resistance
renal_vascular_resistance = preafferent_arteriole_resistance + (afferent_arteriole_resistance + efferent_arteriole_resistance + peritubular_resistance)/baseline_nephrons; 

#################  Renal blood flow
renal_blood_flow_L_min = ((mean_arterial_pressure_MAP - P_venous) / renal_vascular_resistance); 
renal_blood_flow_ml_hr = renal_blood_flow_L_min * 1000 * 60;


#################  Renal Vasculature Pressures
preafferent_pressure = mean_arterial_pressure_MAP - renal_blood_flow_L_min*preafferent_arteriole_resistance;
glomerular_pressure = (mean_arterial_pressure_MAP  - renal_blood_flow_L_min * (preafferent_arteriole_resistance + afferent_arteriole_resistance / baseline_nephrons));
postglomerular_pressure = (mean_arterial_pressure_MAP  - renal_blood_flow_L_min * (preafferent_arteriole_resistance + (afferent_arteriole_resistance+efferent_arteriole_resistance) / baseline_nephrons));


#################  Autoregulatory signals for preafferent and afferent resistances
preaff_autoreg_int = 1 - preaff_autoreg_scale/2;
preafferent_pressure_autoreg_function = preaff_autoreg_int+preaff_autoreg_scale/(1+exp((nom_preafferent_pressure - preafferent_pressure)/myogenic_steepness));
gp_autoreg_int = 1 - gp_autoreg_scale/2;
glomerular_pressure_autoreg_function = gp_autoreg_int+gp_autoreg_scale/(1+exp((nom_glomerular_pressure - glomerular_pressure)/myogenic_steepness));



####################################### Glomerular Filtration #################################################

#################  Glomerular Ultrafiltration Coefficient Kf
glomerular_hydrostatic_conductance_Kf = nom_Kf;
 
#################  Glomerular Filtration Rate
#GFR is calculated according to Starling's equation
number_of_functional_nephrons = baseline_nephrons;
net_filtration_pressure = glomerular_pressure - oncotic_pressure_difference - P_bowmans;
SNGFR_nL_min = glomerular_hydrostatic_conductance_Kf * (glomerular_pressure - oncotic_pressure_difference - P_bowmans);
GFR =  (SNGFR_nL_min / 1000 / 1000000 * number_of_functional_nephrons);
GFR_ml_min = GFR * 1000;

#################  Serum Creatinine
serum_creatinine_concentration = serum_creatinine/blood_volume_L;
creatinine_clearance_rate = GFR_ml_min * dl_ml * serum_creatinine_concentration; #Units: mg/min

#################  Oncotic pressure
#Landis Pappenheimer equation used to calculate oncotic pressure at entrance and exit to glomerulus
#Oncotic pressure is approximated as varying linearly along the glomerulus. Oncotic pressure in the Bowman's space is zero
#Thus the average pressure difference is the average of the entrance and exit oncotic pressure
#We do not consider filtration equilibrium
Oncotic_pressure_in = 1.629*plasma_protein_concentration+0.2935*(plasma_protein_concentration^2);
SNRBF_nl_min = 1e6*1000*renal_blood_flow_L_min/number_of_functional_nephrons;
plasma_protein_concentration_out = SNRBF_nl_min*plasma_protein_concentration/(SNRBF_nl_min-SNGFR_nL_min);
Oncotic_pressure_out = 1.629*plasma_protein_concentration_out+0.2935*(plasma_protein_concentration_out^2);
oncotic_pressure_avg = (Oncotic_pressure_in+Oncotic_pressure_out)/2;


####################################### Plasma sodium concentration and vasopressin secretion #################################################
#################  Plasma sodium concentration
Na_concentration = sodium_amount / blood_volume_L;
ECF_Na_concentration = ECF_sodium_amount/extracellular_fluid_volume;

#################  Control of vasopressin secretion
#A proportional-integral controller is used to ensure there is no steady state error in sodium concentration
#Relative gains of the P and I controller must be chosen carefully.
#In order to permit a steady-state error, the integral controller can be removed. But care should be given then in choosing the proportional gain
Na_water_controller = Na_controller_gain*(Kp_VP*(Na_concentration - ref_Na_concentration)+Ki_VP*Na_concentration_error);

#################  Vasopressin
#Vasopressin is critical in the model, because it allows water excretion to be decoupled from sodium excretion in the collecting duct
normalized_vasopressin_concentration = 1 + Na_water_controller;
vasopressin_concentration = nominal_vasopressin_conc * normalized_vasopressin_concentration;

#Effect of vasopressin on water intake
water_intake_vasopressin_int = 1-water_intake_vasopressin_scale/2;
water_intake = water_intake_species_scale*(nom_water_intake/60/24)*(water_intake_vasopressin_int + water_intake_vasopressin_scale/(1+exp((normalized_vasopressin_concentration_delayed-1)/water_intake_vasopressin_slope)));
daily_water_intake = (water_intake * 24 * 60);


####################################### Tubular Flow and Reabsorption #################################################

Dc_pt = Dc_pt_nom;
L_pt = L_pt_s1_nom+L_pt_s2_nom + L_pt_s3_nom;

#################  Filtered Na Load #################  
SN_filtered_Na_load = (SNGFR_nL_min / 1000 / 1000000)*Na_concentration;
filtered_Na_load = SN_filtered_Na_load*number_of_functional_nephrons;

#################  Regulatory effects on reabsorption #################  

### Tubular Pressure natriuresis effects through RIHP:
pressure_natriuresis_PT_int = 1 - pressure_natriuresis_PT_scale/2;
pressure_natriuresis_PT_effect = max(0.001,pressure_natriuresis_PT_int + pressure_natriuresis_PT_scale / (1 + exp((postglomerular_pressure- RIHP0) / pressure_natriuresis_PT_slope))); 

pressure_natriuresis_LoH_int = 1 - pressure_natriuresis_LoH_scale/2;
pressure_natriuresis_LoH_effect = max(0.001,pressure_natriuresis_LoH_int + pressure_natriuresis_LoH_scale / (1 + exp((postglomerular_pressure - RIHP0) / pressure_natriuresis_LoH_slope))); 

pressure_natriuresis_DCT_magnitude = max(0,pressure_natriuresis_DCT_scale );
pressure_natriuresis_DCT_int = 1 - pressure_natriuresis_DCT_magnitude/2;
pressure_natriuresis_DCT_effect = max(0.001,pressure_natriuresis_DCT_int + pressure_natriuresis_DCT_magnitude/ (1 + exp((postglomerular_pressure - RIHP0) / pressure_natriuresis_DCT_slope))); 

pressure_natriuresis_CD_magnitude = max(0,pressure_natriuresis_CD_scale);
pressure_natriuresis_CD_int = 1 - pressure_natriuresis_CD_magnitude/2;
pressure_natriuresis_CD_effect = max(0.001,pressure_natriuresis_CD_int + pressure_natriuresis_CD_magnitude/ (1 + exp((postglomerular_pressure - RIHP0) / pressure_natriuresis_CD_slope))); 

### AT1-bound AngII effect on PT reabsorption 
AT1_PT_int = 1 - AT1_PT_scale/2;
AT1_effect_on_PT = AT1_PT_int + AT1_PT_scale/(1+exp(-(AT1_bound_AngII - nominal_equilibrium_AT1_bound_AngII)/AT1_PT_slope));

### Aldosterone effect on DCT and CD reabsorption

aldosterone_concentration = normalized_aldosterone_level_delayed* nominal_aldosterone_concentration; 
Aldo_MR_normalised_effect = normalized_aldosterone_level_delayed*MR_antagonist_effect_on_aldo_MR;

aldo_DCT_int = 1 - aldo_DCT_scale/2;
aldo_effect_on_DCT = aldo_DCT_int + aldo_DCT_scale/(1+exp((1 - Aldo_MR_normalised_effect)/aldo_DCT_slope));

aldo_CD_int = 1 - aldo_CD_scale/2;
aldo_effect_on_CD= aldo_CD_int + aldo_CD_scale/(1+exp((1 - Aldo_MR_normalised_effect)/aldo_CD_slope));


### Tubular fractional Reabsorption rates, modulated by regulatory mechanisms (RIHP, AT1-bound AngII, and aldo)
e_pt_sodreab = min(1,nominal_pt_na_reabsorption * AT1_effect_on_PT *pressure_natriuresis_PT_effect);

e_dct_sodreab = min(1,nominal_dt_na_reabsorption * aldo_effect_on_DCT*pressure_natriuresis_DCT_effect *HCTZ_effect_on_DT_Na_reabs); 

e_cd_sodreab = min(1,nominal_cd_na_reabsorption*aldo_effect_on_CD*pressure_natriuresis_CD_effect);

#################  Proximal Tubule Na and Water reabsorption #################  

Na_reabs_per_unit_length = -log(1-e_pt_sodreab)/(L_pt); #mmol/min
Na_pt_out = SN_filtered_Na_load*exp(-Na_reabs_per_unit_length*L_pt);

water_out_pt = ((SNGFR_nL_min / 1000 / 1000000)/SN_filtered_Na_load)*Na_pt_out;

PT_Na_reabs_fraction = 1-Na_pt_out/SN_filtered_Na_load;
PT_water_reabs_fraction = 1-water_out_pt/(SNGFR_nL_min / 1000 / 1000000);

Na_concentration_out_pt = Na_pt_out/water_out_pt;
PT_Na_outflow = Na_pt_out*number_of_functional_nephrons;


#################  Loop of Henle Na and Water reabsorption #################  

##### Descending Loop of Henle #####

water_in_DescLoH = water_out_pt; # L/min
Na_in_DescLoH = Na_pt_out;
Na_concentration_in_DescLoH = Na_concentration_out_pt;

###No solute reabsorption in Descending Limb
Na_out_DescLoH = Na_in_DescLoH;

### Na Reabsorption rate in the Ascending Limb

#The rate of reabsorption per unit length may be flow-dependent, and may be modulated by tubular pressure-natriuresis
# If LoH_flow_dependence = 0, then no flow dependence. If LoH_flow_dependence = 1, perfect flow dependence
deltaLoH_NaFlow = LoH_flow_dependence*(Na_out_DescLoH-nom_Na_in_AscLoH);

# Na reabsorbed per unit length per minute. 
AscLoH_Reab_Rate =(nominal_loh_na_reabsorption*(nom_Na_in_AscLoH+deltaLoH_NaFlow))/L_lh_des; 
effective_AscLoH_Reab_Rate =AscLoH_Reab_Rate*pressure_natriuresis_LoH_effect; 

### Water reabsorption in the Descending Limb
#Descending limb is in equilibrium with interstitium, since water is reabsorbed across the osmotic gradient. The osmotic 
#gradient is generated by Na reabsorption in the ascending limb. Thus, the Na concentration along the descending limb is 
#determined by Na reabsorption in the ascending limb
#Min function necesssary to ensure that the LoH does not reabsorb more Na than is delivered to it
Na_concentration_out_DescLoH = Na_concentration_in_DescLoH*exp(min(effective_AscLoH_Reab_Rate*L_lh_des,Na_in_DescLoH)/(water_in_DescLoH*Na_concentration_in_DescLoH));
water_out_DescLoH = water_in_DescLoH*Na_concentration_in_DescLoH/Na_concentration_out_DescLoH;


##### Ascending Loop of Henle #####
Na_in_AscLoH = Na_out_DescLoH;
Na_concentration_in_AscLoH = Na_concentration_out_DescLoH;
water_in_AscLoH = water_out_DescLoH;

#Na reabsorption along the ascending limb
Na_concentration_out_AscLoH = Na_concentration_in_AscLoH - min(L_lh_des*effective_AscLoH_Reab_Rate, Na_in_DescLoH)*(exp(min(L_lh_des*effective_AscLoH_Reab_Rate, Na_in_DescLoH)/(water_in_DescLoH*Na_concentration_in_DescLoH))/water_in_DescLoH);
Na_reabsorbed_AscLoH = (Na_concentration_in_AscLoH - Na_concentration_out_AscLoH)*water_in_AscLoH;
Na_out_AscLoH = max(0,Na_in_AscLoH - Na_reabsorbed_AscLoH);

#Ascending limb is impermeable to water - no water reabsorption
water_out_AscLoH = water_in_AscLoH;

Na_concentration_out_AscLoH = Na_out_AscLoH/water_out_AscLoH;
LoH_reabs_fraction = 1-Na_out_AscLoH/Na_in_AscLoH;

### Macula Densa Na Flow and Concentration

SN_macula_densa_Na_flow = Na_out_AscLoH;
MD_Na_concentration = Na_concentration_out_AscLoH;

### Tubuloglomerular feedback 
TGF0_tubulo_glomerular_feedback = 1 - S_tubulo_glomerular_feedback/2;
tubulo_glomerular_feedback_signal = (TGF0_tubulo_glomerular_feedback + S_tubulo_glomerular_feedback / (1 + exp((MD_Na_concentration_setpoint - MD_Na_concentration)/ F_md_scale_tubulo_glomerular_feedback)));


#################  Distal Convoluted Tubule #################  

water_in_DCT = water_out_AscLoH; 
Na_in_DCT = Na_out_AscLoH;
Na_concentration_in_DCT = Na_concentration_out_AscLoH; 

#Assume DCT is impermeable to water
water_out_DCT = water_in_DCT;

#Assume sodium reabsorption at a constant fraction of delivery
R_dct = -log(1-e_dct_sodreab)/L_dct;

Na_out_DCT = Na_in_DCT*exp(-R_dct*L_dct);
Na_concentration_out_DCT = Na_out_DCT/water_out_DCT;

DCT_Na_reabs_fraction = 1-Na_out_DCT/Na_in_DCT; 

#################  Collecting Duct #################  

water_in_CD = water_out_DCT;
Na_in_CD = Na_out_DCT;
Na_concentration_in_CD = Na_concentration_out_DCT;

####Assume sodium reabsorbed, then water follows, modulated by vasopressin (ADH)

#Assume sodium reabsorbed at fractional rate eta
R_cd = -log(1-e_cd_sodreab)/L_cd;
Na_out_CD = Na_in_CD*exp(-R_cd*L_cd);
CD_Na_reabs_fraction = 1-Na_out_CD/Na_in_CD; 

#Vasopressin (ADH) effect on water reabsorption through regulation of aquaporin
ADH_water_permeability = min(1,max(0,nom_ADH_water_permeability*normalized_vasopressin_concentration));

#Water reabsorption follows gradient but is regulated by ADH
max_water_reabs_CD = water_in_CD-(Na_concentration_in_CD*water_in_CD-Na_in_CD*(1-exp(-R_cd*L_cd)))/Na_concentration_in_AscLoH; 
water_out_CD = max(0,water_in_CD - ADH_water_permeability*max_water_reabs_CD);

#################  Urine Sodium and Water Excretion #################  

#Urine flow rate
urine_flow_rate = water_out_CD*number_of_functional_nephrons;
daily_urine_flow = (urine_flow_rate * 60 * 24);

#Na Excretion
Na_excretion_via_urine = Na_out_CD*number_of_functional_nephrons;
Na_balance = Na_intake_rate - Na_excretion_via_urine;
water_balance = daily_water_intake - daily_urine_flow;

#Fractional Excretion of Sodium
FENA = Na_excretion_via_urine/filtered_Na_load;


####################################### Tubular Pressure #################################################

#####See written documentation for derivation of the equations below
#flow rates expressed in m3/min, rather than L/min

mmHg_Nperm2_conv = 133.32;
Pc_pt = Pc_pt_mmHg*mmHg_Nperm2_conv;
Pc_lh_des = Pc_lh_des_mmHg*mmHg_Nperm2_conv;
Pc_lh_asc = Pc_lh_asc_mmHg*mmHg_Nperm2_conv;
Pc_dt = Pc_dt_mmHg*mmHg_Nperm2_conv;
Pc_cd = Pc_cd_mmHg*mmHg_Nperm2_conv;
P_interstitial = P_interstitial_mmHg*mmHg_Nperm2_conv;
pi=3.14;

#################  CNT/CD
B1 = (4*tubular_compliance+1)*128*gamma/pi;
mean_cd_water_flow = (water_in_CD-water_out_CD)/2;
B2_cd = (Pc_cd^(4*tubular_compliance))/(Dc_cd^4);
P_in_cd = (0^(4*tubular_compliance+1)+B1*B2_cd*(mean_cd_water_flow/1e3)*L_cd)^(1/(4*tubular_compliance+1));
P_in_cd_mmHg = (P_in_cd+P_interstitial)/mmHg_Nperm2_conv;


#################  DCT
B2_dt = (Pc_dt^(4*tubular_compliance))/(Dc_dt^4);
P_in_dt = (P_in_cd^(4*tubular_compliance+1)+B1*B2_dt*(water_in_DCT/1e3)*L_dct)^(1/(4*tubular_compliance+1));
P_in_dt_mmHg = (P_in_dt+P_interstitial)/mmHg_Nperm2_conv;

#################  Asc LoH
B2_lh_asc = (Pc_lh_asc^(4*tubular_compliance))/(Dc_lh^4);
P_in_lh_asc = (P_in_dt^(4*tubular_compliance+1)+B1*B2_lh_asc*(water_in_AscLoH/1e3)*L_lh_asc)^(1/(4*tubular_compliance+1));
P_in_lh_asc_mmHg = (P_in_lh_asc+P_interstitial)/mmHg_Nperm2_conv;

#################  Desc LoH
A_lh_des = effective_AscLoH_Reab_Rate/(water_in_DescLoH*Na_concentration_in_DescLoH);
B2_lh_des = (Pc_lh_des^(4*tubular_compliance))*(water_in_DescLoH/1e3)/((Dc_lh^4)*A_lh_des);
P_in_lh_des = (P_in_lh_asc^(4*tubular_compliance+1)+B1*B2_lh_des*(1-exp(-A_lh_des*L_lh_des)))^(1/(4*tubular_compliance+1));
P_in_lh_des_mmHg = (P_in_lh_des+P_interstitial)/mmHg_Nperm2_conv;

#################  PT 

A_na = Na_reabs_per_unit_length; 
flow_integral_pt = (SN_filtered_Na_load/A_na)*(1-exp(-A_na*L_pt));

B2_pt = (Pc_pt^(4*tubular_compliance))/(Dc_pt^4);
B3_pt = (SNGFR_nL_min / 1e12)/SN_filtered_Na_load;
P_in_pt= (P_in_lh_des^(4*tubular_compliance+1)+B1*B2_pt*B3_pt*flow_integral_pt)^(1/(4*tubular_compliance+1));
P_in_pt_mmHg = (P_in_pt+P_interstitial)/mmHg_Nperm2_conv;


####################################### Renin Angiotensin Aldosterone System ####################################### 

###Aldosterone is secreted in response to AT1-bound AngII and changes in potassium or sodium concentration
#Potassium concentration is treated as a constant

#AT1-bound AngII effect on Aldosterone
AT1_aldo_int = 1 - AT1_aldo_slope*nominal_equilibrium_AT1_bound_AngII;
AngII_effect_on_aldo = AT1_aldo_int + AT1_aldo_slope*AT1_bound_AngII;

#Normalized aldosterone 
normalized_aldosterone_level = 1*(K_Na_ratio_effect_on_aldo * AngII_effect_on_aldo );

###Renin is secreted in response to decreases in AT1-bound AngII and decreases in MD sodium flow

#Macula Densa Sodium flow effect on renin secretion
#This relationship is known to be non-linear, and md_renin_tau can be calibrated based on data on changes in renin as a functoin of sodium intake
md_effect_on_renin_secretion = md_renin_A*exp(-md_renin_tau*(SN_macula_densa_Na_flow_delayed*baseline_nephrons - nom_LoH_Na_outflow));

#AT1-bound AngII feedback on renin secretion
AT1_bound_AngII_effect_on_PRA = (10 ^ (AT1_PRC_slope * log10(AT1_bound_AngII / nominal_equilibrium_AT1_bound_AngII) + AT1_PRC_yint));

#Aldo effect on renin secretion
aldo_renin_intercept = 1-aldo_renin_slope;
aldo_effect_on_renin_secretion =  aldo_renin_slope * normalized_aldosterone_level + aldo_renin_intercept;

#Plasma renin activity
plasma_renin_activity = concentration_to_renin_activity_conversion_plasma* plasma_renin_concentration*DRI_effect_on_PRA;

#Renin secretion
renin_secretion_rate = (log(2)/renin_half_life)*nominal_equilibrium_PRC*AT1_bound_AngII_effect_on_PRA*md_effect_on_renin_secretion*HCTZ_effect_on_renin_secretion*aldo_effect_on_renin_secretion;

#RAAS degradation rates
renin_degradation_rate = log(2)/renin_half_life; 
AngI_degradation_rate = log(2)/AngI_half_life;
AngII_degradation_rate = log(2)/AngII_half_life;
AT1_bound_AngII_degradation_rate =  log(2)/AT1_bound_AngII_half_life;
AT2_bound_AngII_degradation_rate = log(2)/AT2_bound_AngII_half_life;

#RAAS rate constants
ACE_activity = nominal_ACE_activity*(1 - pct_target_inhibition_ACEi);
chymase_activity = nominal_chymase_activity;
AT1_receptor_binding_rate = nominal_AT1_receptor_binding_rate*(1 - pct_target_inhibition_ARB);
AT2_receptor_binding_rate = nominal_AT2_receptor_binding_rate;


######################################################## ODEs ############################################################# 

#RAAS Pathway
d/dt(AngI) = plasma_renin_activity - (AngI) * (chymase_activity + ACE_activity) - (AngI) * AngI_degradation_rate; 
d/dt(AngII) = AngI * (chymase_activity + ACE_activity) - AngII * AngII_degradation_rate - AngII*AT1_receptor_binding_rate - AngII* (AT2_receptor_binding_rate);
d/dt(AT1_bound_AngII) = AngII * (AT1_receptor_binding_rate) - AT1_bound_AngII_degradation_rate*AT1_bound_AngII;
d/dt(AT2_bound_AngII) = AngII * (AT2_receptor_binding_rate) - AT2_bound_AngII_degradation_rate*AT2_bound_AngII;
d/dt(plasma_renin_concentration) = renin_secretion_rate - plasma_renin_concentration * renin_degradation_rate;

#Change in Extracellular fluid volume over time is determined by the different between water intake and urine outflow
d/dt(blood_volume_L) = C_water_intake_ecf_volume * (water_intake) + C_urine_flow_ecf_volume * (urine_flow_rate) + Q_water*(Na_concentration - ECF_Na_concentration);
d/dt(extracellular_fluid_volume) = Q_water*(ECF_Na_concentration - Na_concentration);

#Change in total body sodium over time is determined by the different between sodium intake and excretion
d/dt(sodium_amount) = C_na_excretion_na_amount * (Na_excretion_via_urine) + C_na_intake_na_amount * (Na_intake_rate) + Q_Na*(ECF_Na_concentration - Na_concentration);
d/dt(ECF_sodium_amount) = Q_Na*(Na_concentration - ECF_Na_concentration);

#These equations serve only to delay the input variable by one timestep. This allows the previous value of the input variable to be used in an equation that appears 
#in the code before the input variable was defined
d/dt(tubulo_glomerular_feedback_effect) = C_tgf * (tubulo_glomerular_feedback_signal-tubulo_glomerular_feedback_effect);
d/dt(normalized_aldosterone_level_delayed) =  C_aldo_secretion * (normalized_aldosterone_level-normalized_aldosterone_level_delayed);
d/dt(preafferent_pressure_autoreg_signal) = 500*(preafferent_pressure_autoreg_function - preafferent_pressure_autoreg_signal);
d/dt(glomerular_pressure_autoreg_signal) = 500*(glomerular_pressure_autoreg_function - glomerular_pressure_autoreg_signal);
d/dt(cardiac_output_delayed) = C_cardiac_output_delayed*(cardiac_output - cardiac_output_delayed);

#Error signals for PI controllers of cardiac output and sodium concentration
d/dt(CO_error) = C_co_error*(cardiac_output-CO_nom);
d/dt(Na_concentration_error) = C_Na_error*(Na_concentration - ref_Na_concentration);

#This equation allows a delay between the secretion of vasopression and its effect on water intake and tubular water reabsorption
d/dt(normalized_vasopressin_concentration_delayed)= C_vasopressin_delay*(normalized_vasopressin_concentration - normalized_vasopressin_concentration_delayed);

#TGF resetting. If C_tgf_reset = 0, no TGF resetting occurs. If it is greater than zero, the setpoint will change over time and will eventually
#come to equal the ambient MD sodium flow rate.
d/dt(F0_TGF) = C_tgf_reset*(SN_macula_densa_Na_flow*baseline_nephrons - F0_TGF);

#As above, these equations allow a variable to be used in equations that appear in the code before the variable was first defined.
d/dt(P_bowmans) = C_P_bowmans*(P_in_pt_mmHg - P_bowmans);
d/dt(oncotic_pressure_difference) = C_P_oncotic*(oncotic_pressure_avg - oncotic_pressure_difference);
d/dt(SN_macula_densa_Na_flow_delayed) = C_md_flow*( SN_macula_densa_Na_flow - SN_macula_densa_Na_flow_delayed);

#Serum Creatinine
d/dt(serum_creatinine) = creatinine_synthesis_rate - creatinine_clearance_rate;

"

save(ode, file = "model_struct.saved")


