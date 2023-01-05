// Description: Bayesian model of retention period assay 
// data : Successive Transfers data; 
//
// infer Pacquire, Pbeta from Successive Transfers data (prob. of acquistion in acq stage and prob of inoculation in inoc stage(s))
// infer Pdeath from Successive Transfers data (assuming same for each transfer stage apart from the initial transfer as qualitatively different)
// infer Pretain from Successive Transfers data (assuming same for each transfer stage)
//  
// note convention throughout is that variables with D_ in front all involve data coming as input

functions {
  vector seq_fun(real start, real end, int N_by)
  {  
    real h;
    vector[N_by]  out;
    h=(end-start) / (N_by-1);
    for (i in 1:N_by)
    { out[i]=start + (i-1)*h; }
    return(out);
  }
}

// DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA 
data {                                    // Data block
int<lower=1> D_Ntransfers;      // number of successive transfers
int<lower=1> D_Wf0;        // Initial number of whitefly in retention groups
int<lower=1> D_RepsR;        // Number of reps for retention experiment
int<lower=0> D_extant1[D_Ntransfers];
int<lower=0> D_infn[D_Ntransfers];                        // Response variable (#if transm)
}

// TRANSFORMED DATA TRANSFORMED DATA TRANSFORMED DATA 
transformed data {                                    // Data block
int<lower=0> D_extant[D_Ntransfers];
for (jj in 1:(D_Ntransfers-1))         // declaring D_extant as want to fit deaths between stages      
D_extant[jj]=D_extant1[jj]-D_extant1[jj+1]; 
D_extant[D_Ntransfers]=D_extant1[D_Ntransfers];
}

// PARAMETERS PARAMETERS PARAMETERS PARAMETERS PARAMETERS 
parameters {                             // Parameters block
real<lower=0,upper=1> compI[1];   
real<lower=0, upper=1> rhoE[1]; 
real<lower=0, upper=1> rhoE0[1]; 
real<lower=0, upper=1> rhoI[1]; 
}


// TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS
transformed parameters {               
  
  real<lower=0,upper=1> chi[D_Ntransfers];
  real<lower=0,upper=1> eta_cum[D_Ntransfers];
  real<lower=0,upper=1> eta[D_Ntransfers];    
  
  chi[1]=rhoE0[1];  // P alive at transfer
  for (jj in 2:(D_Ntransfers))
  chi[jj]=chi[jj-1]*(1-rhoE[1]);
  
  eta[1]=rhoE0[1];  // P alive at transfer
  eta_cum[1]=1-rhoE0[1];  // P die by transfer
  for (jj in 2:(D_Ntransfers))
  {
    eta[jj]=eta[jj-1]*(1-rhoI[1]);
    eta_cum[jj]=eta_cum[jj-1]+rhoE0[1]*((1-rhoI[1])^(jj-2))*rhoI[1];
  }
  
}



// MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL 
model {     
  
  // EXTANT       // EXTANT              // EXTANT              // EXTANT                     
  real matrix_extinct[D_Wf0];  
  
  rhoE[1] ~ beta(1, 1);
  rhoE0[1] ~ beta(1, 1);
  
  for (kk in 1:(D_Ntransfers))
  {
    matrix_extinct = rep_array(0.0,D_Wf0); // note that stage 0 is handled separately in model fitting code
    for (ii in 1:(D_Wf0))
    {
      matrix_extinct[ii]=lchoose(D_Wf0,ii)  // TRANSFER kk+1 extant (depends on rv wf at stage j-1)
      +ii*log(chi[kk])                        /////// log probability extinct @ transfer kk>1 | ii surviving wf
      +(D_Wf0-ii)*log(1-chi[kk])
      +ii*log(rhoE[1]);
    }
    target += D_extant[kk]*log_sum_exp(matrix_extinct);
  }
  
  
  // INFECTION         // INFECTION      // INFECTION      // INFECTION      // INFECTION            
  rhoI[1] ~ beta(1, 1); 
  compI[1] ~ beta(1, 1); 
  
  for (kk in 1:(D_Ntransfers))
  target += binomial_lpmf(D_infn[kk] | D_RepsR, 1-(1-eta[kk]*compI[1])^D_Wf0);
  
}
////////////////////////////////////////////////////////////////////  




// GENERATE GENERATE GENERATE GENERATE GENERATE GENERATE 

generated quantities {      // Generated quantities block. 
real delPeriod[1];
real excessPeriod[1];
real lifespan[1];
real<lower=0, upper=D_RepsR> mnWf[D_Ntransfers]; 
real<lower=0> compI_del[1]; 
real<lower=0, upper=1> pExtant[D_Ntransfers]; 
real<lower=0, upper=D_RepsR> pExtant_Distr[D_Ntransfers]; 
real<lower=0, upper=1> probInfRE[D_Ntransfers]; 
real<lower=0, upper=D_RepsR> probInfSampRE[D_Ntransfers]; 
real<lower=0, upper=1> probCohortExtinctByRE[D_Ntransfers*3]; 

vector[D_Ntransfers] ldeadByVec=seq_fun(0.0,0.0,D_Ntransfers);
vector[3] al_con=[0.25,0.5,0.75]';

delPeriod[1]=(1-rhoI[1])/(1-rhoE[1]);
excessPeriod[1]=rhoI[1]-rhoE[1];
lifespan[1]=1/rhoE[1];

for (jj in 1:(D_Ntransfers))
{
  mnWf[jj]=D_Wf0*chi[jj];
  pExtant[jj]=1-(1-chi[jj])^(D_Wf0);
  pExtant_Distr[jj]=binomial_rng(D_RepsR,pExtant[jj]);
}


for (pepe in 1:3){
  for (tt in 1:(D_Ntransfers)){
    vector[D_Wf0+1] tempyVec=seq_fun(0.0,0.0,D_Wf0+1);
    probInfRE[tt]=1-(1-((1-rhoI[1])^(tt-1))*compI[1])^D_Wf0;
    probInfSampRE[tt]=binomial_rng(D_RepsR,probInfRE[tt]);
    
    // Prob. cohort extinction
    for (kk in 0:(D_Wf0)){
      tempyVec[kk+1]=lchoose(D_Wf0,kk)
      +kk*log(al_con[pepe])
      +(D_Wf0-kk)*log((1-al_con[pepe]))
      +kk*log((eta_cum[tt]));
    }
    ldeadByVec[tt]=log_sum_exp(tempyVec);
    
    probCohortExtinctByRE[tt+(pepe-1)*D_Ntransfers]=1-((exp(ldeadByVec[tt]))^D_RepsR); //binomial_rng(30,probInfRE[tt])
  }   
}

compI_del[1]=compI[1]/rhoE0[1];  

}

