data {

// Robust degree freedom
  real<lower=1> nu;

// Observations
  int<lower=0> N; // Number of observations

// Groups
  vector[N] numG1; // Number of sources. at the moment AZ, clinvar and GWAS-QTL
  vector[N] numG2; // Subgroup indicators qtl type

// Intrecept
  vector[N] bO; // Intercept - Burden test
  vector[N] bOse; // Intercept SE

// // Measure of the X axis - i.e. protein function
// QTLs
  vector[N] xc;
  vector[N] xcse;
// Protein In-silico predictions
//  vector[N] as_conservation;
//  vector[N] as_sift;
//  vector[N] as_polyphen;
  vector[N] as_cadd;
  vector[N] as_alphamissense;
  vector[N] as_revel;
//  vector[N] as_blosum62;
//  vector[N] as_foldx;
//  vector[N] as_consequence;
//  vector[N] as_plddt;

// // Measure of the Y axis - i.e. disease risk

// This are coming from GWAS (including AZ rare-variants one)
  vector[N] yOR; // Response variable
  vector[N] yORse; // Response variable

// Phenotype severity approximation for ClinVar variants
  vector[N] as_clinicalSignificance;
//  vector[N] as_primateai;
}

parameters {
  // Vectors with common variants effects
  vector[N] xcest;
  vector[N] yORest;
  real intercept_random; // Intercept (informed by burden test bO, AZ, and ClinVar rare variants)
  real slope; // Slope for protein function
  vector[5] slope_random; // Random effects for the slope
  vector[N] protein_prior;
  vector[N] disease_prior;
}
transformed parameters {
  // This section if to determine where the information is present - used at the end in the hierarchical model
  real<lower=0> is_eQTL_present = 0;
  real<lower=0> is_pQTL_present = 0;
  real<lower=0> is_CC_present = 0;
  real<lower=0> is_AZ_present = 0;
  real<lower=0> is_CV_present = 0;
  for (n in 1:N) {
    if (numG1[n] == 0 && numG2[n] == 0) {
      is_eQTL_present = 1;
      break;
    }
  }
  for (n in 1:N) {
    if (numG1[n] == 0 && numG2[n] == 1) {
      is_pQTL_present = 1;
      break;
    }
  }
  for (n in 1:N) {
    if (numG1[n] == 3) {
      is_CC_present = 1;
      break;
    }
  }
  for (n in 1:N) {
    if (numG1[n] == 1) {
      is_AZ_present = 1;
      break;
    }
  }
  for (n in 1:N) {
    if (numG1[n] == 2) {
      is_CV_present = 1;
      break;
    }
  }
}
model {
// Protein
// sd have been calculated on the sd of the different tools in the whole prediction set
// protein_prior ~ normal( as_blosum62, .05);
// protein_prior ~ normal( as_foldx, .05);
// protein_prior ~ normal( as_plddt, .05);
// protein_prior ~ normal( as_conservation, .1);
// protein_prior ~ normal( as_sift, .05);
protein_prior ~ normal( as_revel, .28);
protein_prior ~ normal( as_cadd, .13); // for the moment is seems cadd only gives the best outcome - so I may comment the other predictors
protein_prior ~ normal( as_alphamissense, .3);
// protein_prior ~ normal( as_consequence, .1);
// Disease
disease_prior ~ normal( as_clinicalSignificance, .2);
// disease_prior ~ normal( as_primateai, .2); // This is very noisy and doesn't add much to the model

// Priors for latent effect-size parameters — ensures proper posterior
// even when no QTL data is present (e.g. ClinVar-only diseases)
xcest ~ normal(0, 0.2);
yORest ~ normal(0, 1.0);

// Intercept prior — ensures proper posterior when no burden/AZ/ClinVar data present
intercept_random ~ normal(0, 10);
// Slope prior — normal(0, 5) matches single-variant model; slope of ±5 is already extreme (e^5 ≈ 150-fold risk)
slope ~ normal( 0, 5 );
slope_random ~ normal( 0, 5 );
slope_random[5] ~ normal( -1, 5 ); // Rare variants in clinvar are more likely to have negative slope

// Posterior for the intercept
// if not empty use the following bO prior
if (bO[1] != 0) {
  bO ~ normal(intercept_random, bOse);
}

// Measurement models — applied once for all N variants
xc ~ normal( xcest, xcse);
yOR ~ normal( yORest, yORse);

// Posterior for the slope
// Inspiration for this function comes from here: https://dhemery.github.io/DHE-Modules/technical/sigmoid/
// and https://dinodini.wordpress.com/2010/04/05/normalized-tunable-sigmoid-functions/
// The calculation of the sd for the regression comes form Sun et al. 2022 Nature - Genetic associations of protein-coding variants in human disease
for (n in 1:N) {
  if (numG1[n] == 0) { // Common variants
    if ( numG2[n] == 0 ) { // eQTL
      yORest[n] ~ student_t( nu, xcest[n] * slope_random[1], abs(yOR[n]) / fmax(abs(xc[n]), 0.01));
    } else
    if ( numG2[n] == 1 ) { // pQTL
      yORest[n] ~ student_t( nu, xcest[n] * slope_random[2], abs(yOR[n]) / fmax(abs(xc[n]), 0.01));
      }
    }
  else
  if (numG1[n] == 3) { // common coding GWAS
    yORest[n] ~ student_t( nu, protein_prior[n] * slope_random[3], abs(yOR[n] / protein_prior[n]) );
  }
  else
  if (numG1[n] == 1) { // AZ PheWAS
    yORest[n] ~ student_t( nu, intercept_random + slope_random[4] * protein_prior[n], abs(yOR[n] / protein_prior[n]) );
  }
  else
  if (numG1[n] == 2) { // Rare variants
    disease_prior[n] ~ student_t( nu, intercept_random + slope_random[5] * protein_prior[n], 0.3);
  }
}

// Aggregate the slope estimates of the different groups in the hierarchical slope only if the data is present
// If data is not present it will only use the prior - which assume no effect (i.e. no observation = no effect)

if ( is_eQTL_present == 1){
  slope ~ normal( slope_random[1], 0.5 );
}
if ( is_pQTL_present == 1){
  slope ~ normal( slope_random[2], 0.5 );
}
if( is_CC_present == 1){
  slope ~ normal( slope_random[3], 0.5 );
}
if( is_AZ_present == 1){
  slope ~ normal( slope_random[4], 0.5 );
}
if( is_CV_present == 1){
  slope ~ normal( slope_random[5], 0.8 );
}

}
