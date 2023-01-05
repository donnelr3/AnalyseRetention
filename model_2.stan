// Description: Bayesian model of retention period experiment 
// data : Successive Transfers data;
//
// infer Pacquire*Pbeta_j from Successive Transfers data (specific to transfer stage j)
// infer Pdeath from Successive Transfers data (assuming same for each transfer stage apart from the initial transfer as qualitatively different)
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

// DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA 
data {                                    // Data block
int<lower=1> D_Ntransfers;      // number of successive transfers
int<lower=1> D_Wf0;        // Initial number of whitefly in retention groups
int<lower=1> D_RepsR;        // Number of reps for retention experiment
int<lower=0> D_extant1[D_Ntransfers];
int<lower=0> D_infn[D_Ntransfers];                        
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
real<lower=0,upper=1> compI[D_Ntransfers];   // per capita infected insect tranmission probability
real<lower=0, upper=1> rhoE[1]; 
real<lower=0, upper=1> rhoE0[1]; 
}


// TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS 
transformed parameters {               
  
  real<lower=0,upper=1> chi[D_Ntransfers];
  
  chi[1]=rhoE0[1];  // P alive at transfer
  for (jj in 2:(D_Ntransfers))
  chi[jj]=chi[jj-1]*(1-rhoE[1]);
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
  
  for (kk in 1:(D_Ntransfers)){
    compI[kk] ~ beta(1, 1); 
    target += binomial_lpmf(D_infn[kk] | D_RepsR, 1-(1-chi[kk]*compI[kk])^D_Wf0);
  }
  
}
////////////////////////////////////////////////////////////////////  




// GENERATE GENERATE GENERATE GENERATE GENERATE GENERATE 

generated quantities {      // Generated quantities block. 
real<lower=0, upper=D_RepsR> mnWf[D_Ntransfers]; 
real<lower=0> compI_rel[D_Ntransfers];
real<lower=0> compIMaxVal[1];
int<lower=0> loc_compIMaxVal[1];
real<lower=0, upper=1> probExtant[D_Ntransfers]; 
real<lower=0, upper=D_RepsR> probSampExtant[D_Ntransfers]; 
real<lower=0, upper=1> probInfRE[D_Ntransfers]; 
int<lower=0, upper=D_RepsR> probInfSampRE[D_Ntransfers]; 

for (jj in 1:(D_Ntransfers))
{
  mnWf[jj]=D_Wf0*chi[jj];
  probExtant[jj]=1-(1-chi[jj])^(D_Wf0);
  probSampExtant[jj]=binomial_rng(D_RepsR,probExtant[jj]);
}

compIMaxVal[1]=max(compI);
for (jj in 1:(D_Ntransfers)){
  if (compI[jj]==compIMaxVal[1]){
    loc_compIMaxVal[1]=jj;
  }
}

for (tt in 1:(D_Ntransfers)){
  probInfRE[tt]=1-(1-chi[tt]*compI[tt])^D_Wf0;
  probInfSampRE[tt]=binomial_rng(D_RepsR,probInfRE[tt]);
  compI_rel[tt]=compI[tt]/compIMaxVal[1];
}

}

