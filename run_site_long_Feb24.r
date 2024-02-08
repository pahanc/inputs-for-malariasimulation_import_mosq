
Sites_fitted_2021_03_09<-readRDS("Q:\\for_hpc\\Sites_fitted_2021_03_09.RDS")

load("input files\\Admins_in_square.RData")

load("input files\\Level1_names_Af.RData")

load("input files\\Overall_prop_urban.RData")

load("input files\\Rainfall fourier sq params.RData")

load("input files\\Pfpr_sq_2003-18.RData")

load("input files\\Admin_countries.RData")

load("input files\\LLIN IR params.RData")

load("input files\\PBO LLIN IR params.RData")

load("input files\\IRS IR params susc2.RData")

load("input files\\av_IRS_cov_dom_insc.RData")

load("input files\\square_peaks_new.RData")

load("input files\\square_sp_propns_sink.RData")

load("input files\\square bednet distributions_longer.RData")

source("compatibility.R")

run_site<-function(square_number,run_number,total_M,human_population,vaccine_cov,PBO,prop_sp_new,supp_gam,supp_arab,supp_fun){
  
  admins<-which(!is.na(intsct_vec[square_number,]))
  no_admins<-length(admins)
  site_numbers1<-intsct_vec[square_number,admins]
  site_names<-level1_names[site_numbers1]
  site_countries<-admin_country[site_numbers1]
  cat("admin country ", site_countries)
  inds<-which(site_countries=="The Gambia")
  if (length(inds)>0) site_countries[inds]<-"Gambia"
  site_numbers_all<-which(is.element(Sites_fitted_2021_03_09$NAME_1,site_names))
  cat(" database country",Sites_fitted_2021_03_09$NAME_0[site_numbers_all])
  wrong_countries<-which(!is.element(Sites_fitted_2021_03_09$NAME_0[site_numbers_all],site_countries))
  if (length(wrong_countries)>0) site_numbers_all<-site_numbers_all[-wrong_countries]
  
  urban_rural_counts<-table(Sites_fitted_2021_03_09$NAME_1[site_numbers_all])
  if (length(which(urban_rural_counts!=2))>0){# if there are not single urban/rural entries
    #for some sites
    missing_names<-names(urban_rural_counts)[which(urban_rural_counts!=2)]
    site_names_all<-Sites_fitted_2021_03_09[site_numbers_all,"NAME_1"]
    site_info<-cbind(site_names_all,site_numbers_all,Sites_fitted_2021_03_09[site_numbers_all,"ur"])
    colnames(site_info)=c("site_names_all","site_numbers_all","ur_ru")
    k=1
    for (i in 1:length(missing_names)){
      k=which(is.element(site_info[,"site_names_all"],missing_names[i]))
      if (k<nrow(site_info)){
        if (site_info[k,"ur_ru"]=="urban") ur_ru="rural"
        if (site_info[k,"ur_ru"]=="rural") ur_ru="urban"
        vec<-c(site_info[k,"site_names_all"],site_info[k,"site_numbers_all"],ur_ru)
        site_info<-rbind(site_info[1:k,],vec,site_info[(k+1):nrow(site_info),])
      } else {
        site_info<-rbind(site_info,vec)
        }  
    }
    site_numbers_all<-as.numeric(site_info[,"site_numbers_all"])
    site_numbers_urban<-as.numeric(site_info[which(site_info[,"ur_ru"]=="urban"),"site_numbers_all"])
    site_numbers_rural<-as.numeric(site_info[which(site_info[,"ur_ru"]=="rural"),"site_numbers_all"])
  } else {
    site_numbers_urban<-site_numbers_all[which(Sites_fitted_2021_03_09[site_numbers_all,"ur"]=="urban")]
    site_numbers_rural<-site_numbers_all[which(Sites_fitted_2021_03_09[site_numbers_all,"ur"]=="rural")]
  }
  
  params<-get_parameters(square_number=square_number,run_number=run_number,
                         supp_gam=supp_gam,supp_arab=supp_arab,supp_fun=supp_fun,
                         list(prevalence_rendering_min_ages = 2*365,
                              prevalence_rendering_max_ages = 10 * 365,
                              incidence_rendering_min_ages = 0,
                              incidence_rendering_max_ages = 100 * 365,
                              clinical_incidence_rendering_min_ages=c(0,0),
                              clinical_incidence_rendering_max_ages=c(5*365,100*365),
                              severe_incidence_rendering_min_ages= 0*365,
                              severe_incidence_rendering_max_ages= 5*365))
  
  #' Title
  #'
  #' @param params 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  set_site_params<-function(params){
    
    params$model_seasonality = FALSE
    
#***************Find an approximate equilibrium to initialize the model************

    jamie_params_mod<-translate_parameters(params)
    params_mod<-params
    params_mod$species_proportions<-1
    
    #Function to work out the total_M for a give EIR, and the above two parameter sets
    total_M_eir<-function(eir) {
      malariasimulation::set_equilibrium(
        params_mod,
        eir,
        eq_params=jamie_params_mod
      )$total_M[1]
    }
    
    #Find the EIR for the pre-specified total_M
    total_M_eq<-mean(total_M*params$mosq_seasonality[[1]][1:365])
    
    find_total_M<-function(eir){
      return(total_M_eq-total_M_eir(eir))
    }
    
    lower=0.1;upper=600
    while (sign(find_total_M(lower)*find_total_M(upper))>0){
      upper=upper+100
    }
    the_eir<-uniroot(find_total_M,lower=lower,upper=upper,tol=10)
    
    #Set the initial equilibrium state
    params <- malariasimulation::set_equilibrium(params_mod, the_eir$root, eq_params=jamie_params_mod)
    
#***********************Setting vector species params*************************************
#*
    site_species_urban<-sapply(Sites_fitted_2021_03_09[site_numbers_urban,]$vectors,rbind)
    
    site_species_rural<-sapply(Sites_fitted_2021_03_09[site_numbers_rural,]$vectors,rbind)
    
    #Use input data to estimate species proportions
    sp_propns<-square_sp_propns_sink[square_number,c(3,1,2)]
    sp_propns<-c((1-prop_sp_new)*sp_propns,prop_sp_new)
    
    #Vector behaviour parameters
    Q0<-site_species_rural[c("gamb_ss_Q0","arab_Q0","fun_Q0"),1]
    phi_indoors<-site_species_rural[c("gamb_ss_Q_in","arab_Q_in","fun_Q_in"),1]
    phi_bednets<-site_species_rural[c("gamb_ss_Q_bed","arab_Q_bed","fun_Q_bed"),1]
    endophily<-c(0.813,0.422,0.813)
    
    if (length(params$mum)>1) {mum_gamb<-params$mum[1];mum_arab<-params$mum[2];mum_fun<-params$mum[3]
    } else {mum_gamb<-params$mum;mum_arab<-params$mum;mum_fun<-params$mum}
    
    gamb_params<-c('gamb',params$blood_meal_rates,Q0[1],params$foraging_time,phi_bednets[1],phi_indoors[1],mum_gamb)
    names(gamb_params)<-c('species','blood_meal_rates','Q0','foraging_time','phi_bednets','phi_indoors','mum')
    arab_params<-c('arab',params$blood_meal_rates,Q0[2],params$foraging_time,phi_bednets[2],phi_indoors[2],mum_arab)
    names(arab_params)<-c('species','blood_meal_rates','Q0','foraging_time','phi_bednets','phi_indoors','mum') 
    fun_params<-c('fun',params$blood_meal_rates,Q0[3],params$foraging_time,phi_bednets[3],phi_indoors[3],mum_fun)
    names(fun_params)<-c('species','blood_meal_rates','Q0','foraging_time','phi_bednets','phi_indoors','mum') 
    new_params<-c('new',params$blood_meal_rates,Q0[2],params$foraging_time,phi_bednets[2],phi_indoors[2],mum_arab)
    names(new_params)<-c('species','blood_meal_rates','Q0','foraging_time','phi_bednets','phi_indoors','mum') 
    
    params<-set_species(params, species=list(gamb_params,arab_params,fun_params,new_params),proportions=sp_propns)

#***************************Setting human treatment parameters***************************
#*

    burn_in_len=10
    gdrive_len=13

    params <- set_drugs(params, list(AL_params, DHA_PQP_params, SP_AQ_params)) # setting
    
    act_cov_yrs<-Sites_fitted_2021_03_09[site_numbers_all[1],]$interventions[[1]]$year
    
    act_cov_urban<-(matrix(unlist(sapply(Sites_fitted_2021_03_09[site_numbers_urban,]$interventions,cbind)["prop_act",]),ncol=no_admins)* 
                      matrix(unlist(sapply(Sites_fitted_2021_03_09[site_numbers_urban,]$interventions,cbind)["tx",]),ncol=no_admins)) %*% 
      (intsct_pop[square_number,1:no_admins]*urban_prop[site_numbers1])
    
    act_cov_rural<-(matrix(unlist(sapply(Sites_fitted_2021_03_09[site_numbers_rural,]$interventions,cbind)["prop_act",]),ncol=no_admins)* 
                      matrix(unlist(sapply(Sites_fitted_2021_03_09[site_numbers_rural,]$interventions,cbind)["tx",]),ncol=no_admins)) %*% 
      (intsct_pop[square_number,1:no_admins]*(1-urban_prop[site_numbers1]))
    
    act_cov<-act_cov_urban+act_cov_rural
    
    other_drug_cov_urban<-((1-matrix(unlist(sapply(Sites_fitted_2021_03_09[site_numbers_urban,]$interventions,cbind)["prop_act",]),ncol=no_admins))* 
                             matrix(unlist(sapply(Sites_fitted_2021_03_09[site_numbers_urban,]$interventions,cbind)["tx",]),ncol=no_admins)) %*% 
      (intsct_pop[square_number,1:no_admins]*urban_prop[site_numbers1])
    
    other_drug_cov_rural<-((1-matrix(unlist(sapply(Sites_fitted_2021_03_09[site_numbers_urban,]$interventions,cbind)["prop_act",]),ncol=no_admins))* 
                             matrix(unlist(sapply(Sites_fitted_2021_03_09[site_numbers_urban,]$interventions,cbind)["tx",]),ncol=no_admins)) %*% 
      (intsct_pop[square_number,1:no_admins]*(1-urban_prop[site_numbers1]))
    
    other_drug_cov<-other_drug_cov_urban+other_drug_cov_rural
    
    other_drug_cov_yrs<-Sites_fitted_2021_03_09[site_numbers_all[1],]$interventions[[1]]$year
    
    #Add baseline values at the start and repeat end year to 41 years
    act_cov_yrs_warmup<-NULL
    for (i in seq(burn_in_len,1,-1)){
      act_cov_yrs_warmup<-c(act_cov_yrs_warmup,act_cov_yrs[1]-i)
    }
    act_cov_yrs_gdrive<-act_cov_yrs[19]+seq(1,(gdrive_len))
    act_cov_yrs<-c(act_cov_yrs_warmup,act_cov_yrs,act_cov_yrs_gdrive)
    act_cov<-c(rep(0,burn_in_len),act_cov,rep(act_cov[19],(gdrive_len)))
    other_drug_cov_yrs_warmup<-NULL
    for (i in seq(burn_in_len,1,-1)){
      other_drug_cov_yrs_warmup<-c(other_drug_cov_yrs_warmup,other_drug_cov_yrs[1]-i)
    }
    other_drug_cov_yrs_gdrive<-other_drug_cov_yrs[19]+seq(1,(gdrive_len))
    other_drug_cov_yrs<-c(other_drug_cov_yrs_warmup,other_drug_cov_yrs,other_drug_cov_yrs_gdrive)
    other_drug_cov<-c(rep(0,burn_in_len),other_drug_cov,rep(other_drug_cov[19],(gdrive_len)))
    
    params<-set_clinical_treatment(params, 1,((act_cov_yrs-act_cov_yrs[1]+1)*365-365+1),act_cov)

    params<-set_clinical_treatment(params, 3,(other_drug_cov_yrs-other_drug_cov_yrs[1]+1)*365-365+1,other_drug_cov)
    
#***************************Setting SMC parameters***************************
    
    smc_cov_urban<- matrix(unlist(sapply(Sites_fitted_2021_03_09[site_numbers_urban,]$interventions,cbind)["smc",]),ncol=no_admins) %*% 
      (intsct_pop[square_number,1:no_admins]*urban_prop[site_numbers1])
    smc_cov_rural<- matrix(unlist(sapply(Sites_fitted_2021_03_09[site_numbers_rural,]$interventions,cbind)["smc",]),ncol=no_admins) %*% 
      (intsct_pop[square_number,1:no_admins]*(1-urban_prop[site_numbers1]))
    
    smc_cov<-smc_cov_urban+smc_cov_rural
    
    smc_cov_yrs<-Sites_fitted_2021_03_09[site_numbers_all[1],]$interventions[[1]]$year
    
    #Add baseline values at the start
    smc_cov_yrs_warmup<-NULL
    for (i in seq(burn_in_len,1,-1)){
      smc_cov_yrs_warmup<-c(smc_cov_yrs_warmup,smc_cov_yrs[1]-i)
    }
    smc_cov_yrs_gdrive<-smc_cov_yrs[19]+seq(1,(gdrive_len))
    smc_cov_yrs<-c(smc_cov_yrs_warmup,smc_cov_yrs,smc_cov_yrs_gdrive)
    smc_cov<-c(rep(0,burn_in_len),smc_cov,rep(smc_cov[19],(gdrive_len)))
    
    peak <- square_peaks[which(square_peaks[,1]==square_number),2]
    smc_events = data.frame(
      timestep = (smc_cov_yrs-smc_cov_yrs[1]+1) * 365 + peak - 30
    )
    
    params <- set_smc(
      params,
      drug = 3,
      timesteps = smc_events$timestep,
      coverages = smc_cov,
      min_age = 2 * 365 - 1,
      max_age = 5 * 365
    )
    
#***************************Setting bednet parameters***************************
    
    bednet_yrs<-Sites_fitted_2021_03_09[site_numbers_all[1],]$interventions[[1]]$year
    bednet_cov<-distribution_list[[which(squares==square_number)]]
    
    #Add baseline values at the start
    bednet_yrs_warmup<-NULL
    for (i in seq(burn_in_len,1,-1)){
      bednet_yrs_warmup<-c(bednet_yrs_warmup,bednet_yrs[1]-i)
    }
    bednet_yrs_gdrive<-bednet_yrs[19]+seq(1,(gdrive_len))
    bednet_yrs<-c(bednet_yrs_warmup,bednet_yrs,bednet_yrs_gdrive)
    bednet_cov<-c(rep(0,burn_in_len),bednet_cov[1:19],bednet_cov[20:(20+gdrive_len-1)])#,rep(bednet_cov[19],gdrive_len))

    #Option to switch to PBO nets during gene drive releases
    if (!PBO){
     dn04<-dn03;rn04<-rn03;gamman4<-gamman3
    }
    if (PBO){
    dn04<-dn03b;rn04<-rn03b;gamman4<-gamman3b
    }
      
    dn0<-dn03[[square_number]]
    dn0<-c(rep(0.446,15),dn0,rep(dn04[[square_number]][13],gdrive_len+1))
    dn0<-cbind(dn0,dn0,dn0,dn0)

    rn0<-rn03[[square_number]]
    rn0<-c(rep(0.506,15),rn0,rep(rn04[[square_number]][13],gdrive_len+1))
    rn0<-cbind(rn0,rn0,rn0,rn0)

    gamman<-gamman3[[square_number]]*365
    gamman<-c(rep(2.7*365,15),gamman,rep(gamman4[[square_number]][13]*365,gdrive_len+1))
    
    rnm<-rep(0.24,(15+13+gdrive_len+1))
    rnm<-cbind(rnm,rnm,rnm,rnm)

    
    params <- set_bednets(params,timesteps = (bednet_yrs-bednet_yrs[1]+1)*365-365+1, coverages = bednet_cov,retention = 5 * 365, 
                          dn0=as.matrix(dn0), rn=as.matrix(rn0),rnm=as.matrix(rnm),gamman=gamman)
    
#***************************Setting IRS parameters***************************
    irs_yrs<-Sites_fitted_2021_03_09[site_numbers_all[1],]$interventions[[1]]$year
    irs_cov<-av_IRS_list3[[square_number]]
    irs_cov<-c(rep(0,5),irs_cov,irs_cov[length(irs_cov)])
    irs_rounds<-Sites_fitted_2021_03_09[site_numbers_all[1],]$interventions[[1]]$irs_rounds

    #Add baseline values at the start
    irs_yrs_warmup<-NULL
    for (i in seq(burn_in_len,1,-1)){
      irs_yrs_warmup<-c(irs_yrs_warmup,irs_yrs[1]-i)
    }
    irs_yrs_gdrive<-irs_yrs[19]+seq(1,(gdrive_len))
    irs_yrs<-c(irs_yrs_warmup,irs_yrs,irs_yrs_gdrive)
    irs_cov<-c(rep(0,burn_in_len),irs_cov,rep(irs_cov[19],(gdrive_len)))
    
    irs_rounds<-c(rep(0,burn_in_len),irs_rounds,rep(irs_rounds[19],(gdrive_len)))
    
    ls_theta3<-ls_theta3_susc;ls_gamma3<-ls_gamma3_susc;ks_theta3<-ks_theta3_susc;ks_gamma3<-ks_gamma3_susc;ms_theta3<-ms_theta3_susc;ms_gamma3<-ms_gamma3_susc
    
    ls_theta<-ls_theta3[[square_number]]
    ls_theta<-c(rep(ls_theta[1],15),ls_theta,rep(ls_theta[13],(gdrive_len+1)))
    ls_theta<-cbind(ls_theta,ls_theta,ls_theta,ls_theta)
    ls_gamma<-ls_gamma3[[square_number]]
    ls_gamma<-c(rep(ls_gamma[1],15),ls_gamma,rep(ls_gamma[13],(gdrive_len+1)))
    ls_gamma<-cbind(ls_gamma,ls_gamma,ls_gamma,ls_gamma)
    ks_theta<-ks_theta3[[square_number]]
    ks_theta<-c(rep(ks_theta[1],15),ks_theta,rep(ks_theta[13],(gdrive_len+1)))
    ks_theta<-cbind(ks_theta,ks_theta,ks_theta,ks_theta)
    ks_gamma<-ks_gamma3[[square_number]]
    ks_gamma<-c(rep(ks_gamma[1],15),ks_gamma,rep(ks_gamma[13],(gdrive_len+1)))
    ks_gamma<-cbind(ks_gamma,ks_gamma,ks_gamma,ks_gamma)
    ms_theta<-ms_theta3[[square_number]]
    ms_theta<-c(rep(ms_theta[1],15),ms_theta,rep(ms_theta[13],(gdrive_len+1)))
    ms_theta<-cbind(ms_theta,ms_theta,ms_theta,ms_theta)
    ms_gamma<-ms_gamma3[[square_number]]
    ms_gamma<-c(rep(ms_gamma[1],15),ms_gamma,rep(ms_gamma[13],(gdrive_len+1)))
    ms_gamma<-cbind(ms_gamma,ms_gamma,ms_gamma,ms_gamma)

    params <- set_spraying(
      params,
      timesteps = (irs_yrs-irs_yrs[1]+1)*365-365+1,
      coverages = irs_cov, ls_theta = as.matrix(ls_theta),ls_gamma = as.matrix(ls_gamma),
      ks_theta = as.matrix(ks_theta),ks_gamma = as.matrix(ks_gamma),
      ms_theta = as.matrix(ms_theta), ms_gamma = as.matrix(ms_gamma)
    )#this sets sprayingparams$spraying to TRUE
    

  #***************************Setting vaccination parameters***************************
    vacc_cov_yrs_warmup<-NULL
    for (i in seq(burn_in_len,1,-1)){
      vacc_cov_yrs_warmup<-c(vacc_cov_yrs_warmup,2000-i)
    }
    vacc_cov_yrs_gdrive<-2018+seq(1,gdrive_len)
    vacc_cov_yrs<-c(vacc_cov_yrs_warmup,2000:2018,vacc_cov_yrs_gdrive)
    vacc_cov<-c(rep(0,burn_in_len),rep(0,length(2000:2018)),rep(vaccine_cov,gdrive_len))
    
    month1<-30;year1<-365
    
    vacc_events = data.frame(
      timestep = (vacc_cov_yrs-vacc_cov_yrs[1]+1) * 365 + peak - 30
    )

     params <- set_rtss_epi(
       params,
       start = vacc_events$timestep[30],
       end = vacc_events$timestep[length(vacc_events$timestep)],
       coverage =vaccine_cov,
       min_wait = 6 * month1,
       age = 5 * month1,
       boosters = c(365),#One booster dose before transmission season
       booster_coverage = c(0.8),
       seasonal_boosters = FALSE
     )

    
  #***********************Set up human population************************************
    
    average_age=22*365 
    
    #Set human population size
    params$individual_mosquitoes=FALSE
    params$human_population<-human_population

    return(params)
  }
  
  params<-set_site_params(params)
  params$total_M<-total_M
  params$human_population<-human_population
  
  output<-run_simulation(365*42,params)
  
  print(paste("species_proportions",params$species_proportions))

  return(list(output,params$species_proportions))
}
