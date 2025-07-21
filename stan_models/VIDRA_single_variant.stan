data {
// Hyperparameters
  real<lower=0> h1;
// Observations
  int<lower=0> N; // Number of observations
  // int<lower=0> K; // Number of total G1 (i.e. number of sources)  
// Groups
  real numG1; // Number of sources. at the moment AZ, clinvar and GWAS-QTL
  real numG2; // Subgroup indicators qtl type
// // Measure X - protein function
// QTLs
  real xc; 
  real xcse; 
// Protein predictions
  real as_conservation;
  real as_sift;
  real as_polyphen;
  real as_cadd;
  real as_alphamissense;
// // Measure Y - disease risk
// This are coming from GWAS (including AZ rare-variants one)
  real yOR; // Response variable 
  real yORse; // Response variable 
// Phenotype predictions
  real as_clinicalSignificance;
  real as_primateai;
}
parameters {
  // Vectors with common variants effects
  real xcest;
  real yORest; 
  real slope; // Slope for protein function
  real protein_prior;
  real disease_prior; 
}
model {
// Protein
protein_prior ~ normal( as_conservation, .05); 
protein_prior ~ normal( as_cadd, .05); 
protein_prior ~ normal( as_alphamissense, .05);
// Disease
disease_prior ~ normal( as_clinicalSignificance, .05);
disease_prior ~ normal( as_primateai, .05);
slope ~ normal( 0, 5 );
if (numG1 == 0) { // Common variants
  // // Priors
  xc ~ normal( xcest, xcse);
  yOR ~ normal( yORest, yORse);
  slope ~ normal(yOR / xc, abs(yORse/xcest));
  } 
else 
if (numG1 == 1) { // common coding GWAS
  slope ~ normal(yORest / protein_prior, abs(yORse/0.1));
} 
else 
if (numG1 == 2) { // AZ PheWAS
  slope ~ normal(yORest / protein_prior, abs(yORse/0.1));
} 
else 
if (numG1 == 3) { // Rare variants
  slope ~ normal(disease_prior / protein_prior, 0.1); 
}
}

