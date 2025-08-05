import sys
# import os
import pandas as pd
import json
import requests
from functools import reduce
import multiprocessing  
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype
from itertools import repeat
from statistics import median
import resource
# import dask.dataframe as dd
# from math import ceil
import asyncio
import aiohttp
import numpy as np
from sys import argv
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import pyarrow as pa
import pyarrow.parquet as pq
import httpstan
import httpstan.models as hst
import stan
# from cmdstanpy import CmdStanModel
# import cmdstanpy
from stan.model import Model
import multiprocessing
from math import floor


TIME_PERIOD = 3600 # 1 hour
REQUESTS_PER_PERIOD = 55000  # Number of requests I can place in 1 hour
MAX_CONCURRENT_REQUESTS = 1000  # Set the maximum number of concurrent requests
# convert request from DF to suitable for ProtVarAPI
def structure_request_ProtVar( dat_frame, chrom, pos, ref, alt ):
    c = [ str(i) for i in dat_frame[chrom] ]
    p = [ str(i) for i in dat_frame[pos] ]
    r = [ str(i) for i in dat_frame[ref] ]
    a = [ str(i) for i in dat_frame[alt] ]
    # zip  - aggregate the df to list of rows
    zip_obj = zip( c, p, r, a )
    # return
    return( [" ".join(i) for i in list(zip_obj)] )

# Reshape the request of ProtVar API and return dummy structured df
def structure_features(pdDF_in):
    try:
        df_out = ( pd.get_dummies( 
                pdDF_in.reset_index(drop=True).dropna(),
                prefix="",
                prefix_sep='',
                dtype=int 
                ).apply(
                    func = sum, 
                    axis=1, 
                    result_type='broadcast'
                    ).iloc[[0]] 
            )
        return(df_out)
    except:
        pass

async def fetch_protvar_data(session, chrom, pos, ref, alt, semaphore):
    async with semaphore:
        variant = f"{chrom} {pos} {ref} {alt}"
        headers = {
            'accept': 'application/json',
            'Content-Type': 'application/json',
        }
        params = {
            'function': 'true',
            'population': 'true',
            'structure': 'true',
        }
        url = 'https://www.ebi.ac.uk/ProtVar/api/mappings'      
        attempts = 0
        while attempts < 10:
            try:
                async with session.post(url, params=params, headers=headers, json=[variant], timeout=20) as response:
                    return await response.json()
            except aiohttp.ClientError as e:
                attempts += 1
                print(f"Warning: {e}. Re-entering the loop for now ({10 - attempts} attempts left)")
            except asyncio.TimeoutError as e2:
                attempts += 1
                print(f"Warning: {e2}. time-out error. Re-entering the loop for now ({10 - attempts} attempts left)")
        return None

async def process_PVvariant(chrom, pos, ref, alt, semaphore):
    async with aiohttp.ClientSession() as session:
        response_data = await fetch_protvar_data(session, chrom, pos, ref, alt, semaphore)   
    if response_data is None or len(response_data) == 0:
        print("response_data EMPTY")
        return pd.DataFrame()
    # Process the response_data and extract relevant information
    # Additional filter if variant return invalid entry
    try:
        if (pd.json_normalize(response_data['inputs'][0]['messages']).type[0]) == 'ERROR':
            print(f"Warning: {chrom}_{pos}_{ref}_{alt} returned error message from ProtVar")
            return( pd.DataFrame() )
    except:
        pass
    
    try:
        pd.json_normalize(response_data['inputs'])
    except KeyError as e:
        print ('I got a KeyError - reason "%s"' % str(e))
        return( pd.DataFrame() )
    
    # Additional filter if variant return empty entry
    empty_r_check = pd.json_normalize( response_data, record_path=['inputs','mappings','genes'])
    if ( empty_r_check.empty ):
        print(f"Warning: {chrom}_{pos}_{ref}_{alt} returned empty response from ProtVar")
        return( pd.DataFrame() )
    
    # Give a first parse of the output
    # Unfortunately the output is not standard across all the entries - meaning that some would need manual curation
    general_info = pd.json_normalize(response_data['inputs'])
    gene_info = pd.json_normalize(
        response_data, 
        record_path= ['inputs','mappings','genes'], 
        errors='ignore')
    isoform_info = pd.json_normalize(
        response_data, 
        record_path= ['inputs','mappings','genes','isoforms'], 
        errors='ignore')
    iso_canonical = isoform_info.loc[isoform_info.canonical]
    
    description = pd.DataFrame()
    type_feature = pd.DataFrame()
    category = pd.DataFrame()
    typeDescription = pd.DataFrame()
    
    try:
        for i in iso_canonical['referenceFunction.features']:
            tmp = pd.json_normalize(i)
            description = pd.concat([description, tmp.description],axis=0)
            type_feature = pd.concat([type_feature, tmp.type],axis=0)
            category = pd.concat([category, tmp.category],axis=0)
            typeDescription = pd.concat([typeDescription, tmp.typeDescription],axis=0)
    except:
        description = pd.DataFrame()
        type_feature = pd.DataFrame()
        category = pd.DataFrame()
        typeDescription = pd.DataFrame()
    
    feat_desctiption = structure_features(description)
    feat_type = structure_features(type_feature)
    feat_category = structure_features(category)
    feat_typeDescription = structure_features(typeDescription)
    
    try:
        iso_canonical['referenceFunction.comments']
    except:
        return( pd.DataFrame() )
    
    location_evidences = pd.DataFrame()
    for i in iso_canonical['referenceFunction.comments']:
        tmp_entry_comment = pd.json_normalize(i)
        for _,entry in tmp_entry_comment.iterrows():
            try:
                tmp = pd.json_normalize( entry.locations )
                location_evidences = pd.concat([location_evidences,tmp])
            except:
                pass
    
    try:
        feat_location_value = structure_features(location_evidences['location.value'])
    except:
        feat_location_value = pd.DataFrame()
    
    try:
        feat_location_topol = structure_features(location_evidences['topology'])
    except:
        feat_location_topol = pd.DataFrame()
    
    interaction_evidences = pd.DataFrame()
    for i in iso_canonical['referenceFunction.comments']:
        tmp_entry_comment = pd.json_normalize(i)
        for _,entry in tmp_entry_comment.iterrows():
            try:
                tmp = pd.json_normalize( entry.interactions )
                interaction_evidences = pd.concat([interaction_evidences,tmp])
            except:
                pass
    
    try:
        n_interactors = pd.DataFrame(data=[len(interaction_evidences.accession2.drop_duplicates())], columns=["n_interactors"],index=[0])
    except:
        n_interactors = pd.DataFrame()
    
    function_type = list()
    dbReferences = list()
    for i in iso_canonical['referenceFunction.comments']:
        tmp_entry_comment = pd.json_normalize(i)
        for _,entry in tmp_entry_comment.iterrows():
            try:
                if pd.notna(entry.type):
                    function_type.append( entry.type )
            except:
                pass
            try:
                if pd.notna(entry['reaction.dbReferences']):
                    dbReferences.append( entry['reaction.dbReferences'] )
            except:
                pass
            
    function_typeDF = pd.DataFrame(columns=[*set(function_type)]) 
    try: 
        function_typeDF.loc[0] = 1 
    except: 
        pass
    
    dbReferencesDF = pd.DataFrame(columns=[*set(dbReferences)]) 
    try: 
        dbReferencesDF.loc[0] = 1 
    except: 
        pass
    
    ref_dbReferences_info = pd.DataFrame()
    for i in iso_canonical['referenceFunction.dbReferences']:
        tmp = pd.json_normalize(i)
        ref_dbReferences_info = pd.concat([ref_dbReferences_info,tmp])
    
    try:
        feat_dbReferences_info = structure_features(ref_dbReferences_info.id)
        n_fam_domains = pd.DataFrame(data=[len(ref_dbReferences_info.id)], columns=["n_families"],index=[0])
    except:
        feat_dbReferences_info = pd.DataFrame()
        n_fam_domains = pd.DataFrame()
    
    ref_pockets_info = pd.DataFrame()
    for i in iso_canonical['referenceFunction.pockets']:
        tmp = pd.json_normalize(i)
        ref_pockets_info = pd.concat([ref_pockets_info,tmp])
        
    ref_foldxs_info = pd.DataFrame()
    for i in iso_canonical['referenceFunction.foldxs']:
        tmp = pd.json_normalize(i)
        ref_foldxs_info = pd.concat([ref_foldxs_info,tmp])
    
    ref_interactions_info = pd.DataFrame()
    for i in iso_canonical['referenceFunction.interactions']:
        tmp = pd.json_normalize(i)
        ref_interactions_info = pd.concat([ref_interactions_info,tmp])
    
    pdock_interactions = pd.DataFrame()
    try:
        pdock_interactions = pd.DataFrame(data=[median(ref_interactions_info.pdockq)], columns=["pdock_median"],index=[0])
    except:
        pass
    
    # Get all info for variants that possibly have multiple entries in these features
    populationFrequencies = pd.DataFrame()
    predictions = pd.DataFrame()
    clinicalSignificances = pd.DataFrame()
    for i in iso_canonical['populationObservations.proteinColocatedVariant']:
        tmp = pd.json_normalize(i)
        try:
            tmp = tmp[ tmp['genomicLocation'].str.contains( str(general_info.pos[0])+str(general_info.ref[0])+">"+str(general_info.alt[0]) ) ]
        except:
            tmp = pd.json_normalize(i)
        try:
            for j in tmp.populationFrequencies:
                tmp_freq = pd.json_normalize(j)
                populationFrequencies = pd.concat([populationFrequencies,tmp_freq])
        except:
            pass
        # Becasue I'm looking at aa modification, the same aa variant can be caused by multiple nt variants
        # This lead to multiple lines and multiple redictions here.
        # Hence need to filter the line for the variants I did use as input
        try:
            for j in tmp.predictions:
                tmp_pred = pd.json_normalize(j[0])
                predictions = pd.concat([predictions,tmp_pred])
        except:
            pass
        try:
            for j in tmp.clinicalSignificances:
                tmp_sign = pd.json_normalize(j[0])
                clinicalSignificances = pd.concat([clinicalSignificances,tmp_sign])
        except:
            pass
        
    populationFrequenciesDF = pd.DataFrame()
    try:
        populationFrequenciesDF = pd.DataFrame({ 'MAF': populationFrequencies.frequency }, index=[0])
    except:
        pass
    
    predictionsDF = pd.DataFrame()
    try:
        predictionsDF = pd.DataFrame({
            'prediction_tool': predictions.predAlgorithmNameType, 
            'predScore': predictions.score }).T
        predictionsDF.columns = predictionsDF.iloc[0]
        predictionsDF = predictionsDF.drop(predictionsDF.index[0])
        predictionsDF = predictionsDF.reset_index(drop=True)
    except:
        pass
    
    clinicalSignificancesDF = pd.DataFrame()
    try:
        clinicalSignificancesDF = pd.DataFrame({ 'clinicalSignificances': clinicalSignificances.type }, index=[0]).pipe(pd.get_dummies,dtype='int')
    except:
        pass
    
    protvar_info = pd.concat(
        [ general_info.drop(['messages','mappings','errors','groupBy','type','genType','id'],axis=1),
            gene_info.drop(['isoforms'],axis=1),
            iso_canonical.drop(
                    ['translatedSequences','proteinStructure','populationObservations.genomicColocatedVariant',
                    'populationObservations.proteinColocatedVariant','referenceFunction.geneNames',
                    'referenceFunction.features','referenceFunction.comments','referenceFunction.sequence.sequence',
                    'referenceFunction.sequence.modified','referenceFunction.name','referenceFunction.alternativeNames',
                    'referenceFunction.lastUpdated','referenceFunction.dbReferences','referenceFunction.id',
                    'referenceFunction.proteinExistence','referenceFunction.type',
                    'referenceFunction.position','referenceFunction.accession',
                    'referenceFunction.pockets','referenceFunction.foldxs','referenceFunction.interactions'],
                axis=1, errors='ignore'),
            feat_desctiption,
            feat_type,
            feat_category,
            feat_typeDescription,
            populationFrequenciesDF,
            pdock_interactions,
            predictionsDF, 
            ref_pockets_info.drop(['structId','residList'],axis=1, errors='ignore'),
            ref_foldxs_info,
            dbReferencesDF,
            n_fam_domains,
            feat_location_value,
            feat_location_topol,
            feat_dbReferences_info,
            n_interactors,
            clinicalSignificancesDF,
            function_typeDF
        ], axis=1 )
    # Remove possile duplicated columns
    protvar_info = protvar_info.loc[:,~protvar_info.columns.duplicated()].copy()    
    return(protvar_info)

async def PV_vars(variants):  
    semaphore = asyncio.Semaphore(MAX_CONCURRENT_REQUESTS)  # Control the number of concurrent requests
    tasks = []
    for variant in variants:
        tasks.append(process_PVvariant(variant["chrom"], variant["pos"], variant["ref"], variant["alt"], semaphore))
    results = await asyncio.gather(*tasks)
    final_dataframe = pd.concat(results, ignore_index=True)
    return(final_dataframe)

async def fetch_vep_data(session, chrom, pos, ref, alt, input_type, semaphore):
    async with semaphore:
        if input_type == "gen_coords":
            var = f"{chrom}:g.{pos}{ref}>{alt}"
        elif input_type == "HgvsId":
            var = chrom
        headers = {'Content-type': 'application/json'}
        url = (
            f'https://rest.ensembl.org/vep/human/hgvs/{var}'
            '?EVE=1&Blosum62=1&tsl=1&numbers=1&mutfunc=1&ga4gh_vrs=1'
            '&domains=1&dbscSNV=1&dbNSFP=1&canonical=1&UTRAnnotator=1'
            '&SpliceAI=1&NMD=1&MaveDB=1&LoF=1&IntAct=1&Conservation=1&CADD=1'
            '&DisGeNET=1&Enformer=1&LoF=1&MaveDB=1'
        )
        attempts = 0
        while attempts < 10:
            try:
                async with session.get(url, headers=headers, timeout=30) as response:
                    return await response.json()
            except aiohttp.ClientError as e:
                attempts += 1
                print(f"Warning: {e}. Re-entering the loop for now ({10 - attempts} attempts left)")
            except asyncio.TimeoutError as e2:
                attempts += 1
                print(f"Warning: {e2}. time-out error. Re-entering the loop for now ({10 - attempts} attempts left)")
        return None

async def process_VEPvariant(chrom, pos, ref, alt, input_type, semaphore):
    async with aiohttp.ClientSession() as session:
        response_data = await fetch_vep_data(session, chrom, pos, ref, alt, input_type, semaphore)   
    if response_data is None or not isinstance(response_data, list) or len(response_data) == 0:
        return pd.DataFrame()
    # Extract and process data from response_data
    # VEP Info
    # General VEP info
    vep = ( 
        pd.json_normalize(response_data)[[
            'seq_region_name',
            'allele_string',
            'start',
            'strand',
            'most_severe_consequence',
            'input',
            'assembly_name',
            'end']].set_index('input')
        )
    # Info on colocated variants
    try:
        colocated_variants = pd.json_normalize(
            response_data,  
            record_path=["colocated_variants"],
            meta='input',
            errors='ignore').set_index('input')
        # this bit counts the null entries in the rows and keeps only the most complete row
        colocated_variants['count'] = pd.isnull(colocated_variants).sum(1)
        colocated_variants=colocated_variants.sort_values(['count']).drop(['count'],axis=1).iloc[[0]]
    except:
        colocated_variants = pd.DataFrame()
    # Info on transcript_consequences variants
    try:
        transcript_consequences = (
            pd.json_normalize(
            response_data,
            record_path=["transcript_consequences"],
            meta='input',
            errors='ignore')
            # Get only info on canonical transcript
            .loc[ lambda x: x.canonical == 1 ]
            .set_index('input')
            .explode("consequence_terms")
            )
        # this bit counts the null entries in the rows and keeps only the most complete row
        transcript_consequences['count'] = pd.isnull(transcript_consequences).sum(1)
        transcript_consequences=transcript_consequences.sort_values(['count']).drop(['count'],axis=1).iloc[[0]] 
    except:
        transcript_consequences=pd.DataFrame()
    # These are a subset of transcript_consequences variants
    # extract in a second time because it is easier
    try:
        intact_info = (
            pd.json_normalize(
            transcript_consequences.intact[0],
            errors='ignore')
            .feature_type
            .pipe(pd.get_dummies,dtype="int")
            .agg( sum, result_type='broadcast' )
            .iloc[[0]]
            )
    except:
        intact_info = pd.DataFrame()
    # Extract domain info
    # Same steps of previous request
    try:
        domain_info = (
            pd.json_normalize(
            transcript_consequences.reset_index().domains[0],
            errors='ignore'
            ) 
            .pipe(lambda df_,x : df_[~df_.db.str.startswith(x)], ("ENSP","Gene3D","Alpha") )
            .name
            .pipe(pd.get_dummies,dtype="int")
            .apply( sum, axis=1,result_type='broadcast')
            .iloc[[0]]
            )
    except:
        domain_info = pd.DataFrame()
    # Put everything together - This should be one row per variants
    # Drip intact column if exist - This is because it is better parsed in the intact info df
    transcript_consequences = transcript_consequences.drop(['intact','domains'],axis=1, errors='ignore')
    df = pd.concat(  
        [vep.reset_index(drop=True), 
            colocated_variants.reset_index(drop=True), 
            transcript_consequences.reset_index(drop=True),
            intact_info.reset_index(drop=True),
            domain_info.reset_index(drop=True)
            ], axis=1 )
    # Remove duplicated columns
    df = df.loc[:,~df.columns.duplicated()].copy()
    # attach a variants ID that is good to aggregate with ProtVar
    df['inputStr'] = f"{chrom} {pos} {ref} {alt}"    
    # return output data frame
    return(df)

## Here I can implement a function to ger iteration done out of missing. 
# I can use tqdm to get a progress bar
async def VEP_vars(variants):
    semaphore = asyncio.Semaphore(MAX_CONCURRENT_REQUESTS)
    tasks = []
    for variant in variants:
        tasks.append(process_VEPvariant(variant["chrom"], variant["pos"], variant["ref"], variant["alt"], variant["input_type"], semaphore))
    results = await asyncio.gather(*tasks)
    final_dataframe = pd.concat(results, ignore_index=True)
    return(final_dataframe)

##############################################################################################
##############################################################################################
##############################################################################################
# Function for the estimate of the VIDRA effect

def transform_foldx_score(foldx_score):
    max_score = np.max(foldx_score)
    inverted_score = max_score - foldx_score
    max_inverted_score = np.max(inverted_score)
    transformed_score = 1 - (inverted_score / max_inverted_score)
    return transformed_score

# This function scales the Cadd Score as linear. 
def caddScaled(cadd_score):
    # logarithmically scaled Cadd score typically ranges from 1 to over 100
    linear_score = cadd_score
    min_score, max_Score = 0, 50 # 50 is arbitrary set on my experience
    scaled_score = (linear_score - min_score) / (max_Score - min_score)
    return scaled_score

### Shape Stan vectors
def formatData_StanModel(df, N_Yfeat, N_Xfeat):
    dicDt = {
        # Number of observations
        'N': df.shape[0],
        # Number of feat for Y
        'N_Yfeat': N_Yfeat,
        # Number of feat for X
        'N_Xfeat': N_Xfeat,
        # GROUPS
        #########################
        # Upper bound genetic sources
        'numGsource': len(df.var_type.unique()), 
        # Group indicators
        # to keep always the same order I set categories manually
        # + 1 is to compensate python 0 index
        'GsourceLab': (pd.Categorical( df.var_type, categories=['common_NC', 'common_C', 'rare']).factorize()[0] + 1).tolist(),
        # 'GsourceLabIndex': pd.Categorical(df.var_type,categories=['common_NC', 'common_C', 'rare']).factorize()[1].to_list(),
        #########################
        # Upper bound bio_feature groups
        'numGcell': len(df.bio_feature.unique()), 
        # Group indicators
        # + 1 is to compensate python 0 index
        'GcellLab': (df.bio_feature.factorize()[0] + 1).tolist(),
        #########################
        # Upper bound QTL groups
        'numGqtl': len(df.type_id.unique()),
        # Group indicators
        # + 1 is to compensate python 0 index
        'GqtlLab': (pd.Categorical(df.type_id.fillna(0), categories=[0, 'sqtl', 'eqtl', 'pqtl']).factorize(use_na_sentinel = False)[0] + 1).tolist(),
        # 'GqtlLabIndex': pd.Categorical(df.type_id.fillna(0), categories=[0, 'sqtl', 'eqtl', 'pqtl']).factorize(use_na_sentinel = False)[1].to_list(),
        #########################
        # Upper bound allelic requirement groups
        'numGas': len(df.allelicRequirements.unique()),
        # Group indicators
        # + 1 is to compensate python 0 index
        'GasLab': ( pd.Categorical( df.allelicRequirements.fillna(0), categories=["{'list': array([{'element': 'dominant'}], dtype=object)}", "{'list': array([{'element': 'recessive'}], dtype=object)}"]).factorize(use_na_sentinel = False)[0] + 1).tolist(),
        #########################
        # Upper bound allelic requirement groups
        'numGmod': len(df.Model.unique()),
        # Group indicators
        # + 1 is to compensate python 0 index
        'GmodLab': ( pd.Categorical( df.Model.fillna(0), categories=[0,'allelic', 'dominant', 'recessive']).factorize()[0] + 1).tolist(),
        # Values
        # GWAS - QTL - i.e. common
        # Y variable common
        'yc': df.beta_gwas.fillna(0).to_list(),
        # Y se variable common
        'ycse': df.se_gwas.fillna(0).to_list(),
        # X variable common
        'xc': df.beta_qtl.fillna(1).to_list(),
        # X se variable common
        'xcse': df.se_qtl.fillna(1).to_list(),
        # AZ 
        # Y variable 
        # The filling aggregates the wo columns
        'yOR': df.Oddsratio.fillna( np.exp(df.beta_gwas) ).fillna(df.se_gwas).fillna(1).to_list(),
        'yORse': ((df.OddsratioUCI - df.OddsratioLCI) / 3.92).fillna(df.se_gwas).fillna(1).to_list(), 
        #########################
        # Intercept information - 
        'bO' : df.oddsRatio.fillna(0).to_list(),
        'bOse': ( np.sqrt(df.studyCases) * (df.oddsRatioConfidenceIntervalUpper - df.oddsRatioConfidenceIntervalLower)/3.92).fillna(0).to_list(),
        ##########################
        #### Protein function - x axis
        # normalise between 0 and 1 - where 0 is detremental and 1 if stabilising the variants
        'as_blosum62' : (1 / (1 + np.exp(-df.blosum62))).fillna(1).to_list(),
        'as_foldx' : transform_foldx_score(df.foldxDdq).fillna(1).to_list(), 
        'as_plddt': (df.plddt / 100).fillna(0).to_list(), # scale from 0 to 1 
        'as_conservation' : (1 - df['referenceFunction.conservScore']).fillna(1).to_list(),
        'as_sift' : df.sift_score.fillna(1).to_list(),
        'as_polyphen' : df.polyphen_score.fillna(1).to_list(),
        ## Categories for variants
        'as_consequence': pd.Categorical(
        df.most_severe_consequence.fillna(0),
        categories=[0, 
                '3_prime_UTR_variant', 
                '5_prime_UTR_variant', 
                'splice_region_variant', 
                'splice_donor_region_variant',
                'non_coding_transcript_exon_variant', 
                'synonymous_variant', 
                'splice_donor_variant',
                'missense_variant',
                'stop_gained'])\
        .factorize()[0].tolist(), 
        'as_consequenceIndex': pd.Categorical(
        df.most_severe_consequence.fillna(0),
                categories=[0, 
                        '3_prime_UTR_variant', 
                        '5_prime_UTR_variant', 
                        'splice_region_variant', 
                        'splice_donor_region_variant',
                        'non_coding_transcript_exon_variant', 
                        'synonymous_variant', 
                        'splice_donor_variant',
                        'missense_variant',
                        'stop_gained'])\
                .factorize()[1].to_list(),
        ##########################
        #### Diseases - y axis
        'as_clinicalSignificance': pd.Categorical(
        df.clinicalSignificances.fillna(0),
        categories=[0, 
                'benign', 
                'likely pathogenic', 
                'pathogenic'])\
        .factorize()[0].tolist(),
        'as_clinicalSignificanceIndex': pd.Categorical(
        df.clinicalSignificances.fillna(0),
        categories=[0, 
                'benign', 
                'likely pathogenic', 
                'pathogenic'])\
        .factorize()[1].to_list(),
        'as_cadd': caddScaled(df.caddScore.fillna(0)).to_list(),
        'as_eve': df.eve_score.fillna(0).to_list()
    }
    return(dicDt)


def CleanVEPoutput(vep_file):
    listDF = []
    with open(vep_file, 'r') as file:
        for line in file:
            munchObj = Munch.fromJSON(line)
            # standard info for the variants - i.e. commonc across all transcripts
            if munchObj.variant_class != "SNV":
                continue
            variantInfo = pd.DataFrame({
                'variantId': munchObj.input.replace(". ","").replace(" ","_"),
                'most_severe_consequence': munchObj.most_severe_consequence
            }, index=[0])
            transcriptConsequencesDF = pd.json_normalize(munchObj.transcript_consequences) # This is all parsed correctly in this way
            if transcriptConsequencesDF.columns.str.contains('canonical').any():
                pass
            else:
                print( munchObj.input, 'doesn t have canonical transcript')
                continue
            transcriptConsequencesDF = transcriptConsequencesDF[transcriptConsequencesDF.canonical == 1]
            # get the line with hte list missing information
            transcriptConsequencesDF['missingInfo'] = transcriptConsequencesDF.apply(lambda x: x.isna().sum(), axis=1)
            transcriptConsequencesDF = transcriptConsequencesDF[ transcriptConsequencesDF.missingInfo == transcriptConsequencesDF.missingInfo.min()]
            colKeep = ['primateai', 'revel', 'cadd_phred', 'cadd_raw', 'pli_gene_value', 
                       'impact', 'gene_symbol','conservation', 'loftool', 'canonical', 
                       'consequence_terms', 'polyphen_prediction', 'am_pathogenicity', 
                       'sift_prediction', 'sift_score', 'am_class', 'blosum62', 
                       'polyphen_score']
            transcriptConsequencesDF = transcriptConsequencesDF.filter(colKeep).reset_index(drop=True)
            # aggregate all the information
            tmpDF = pd.concat( [variantInfo, transcriptConsequencesDF], axis=1)
            if tmpDF.shape[0] > 1:
                print('multiple eligible variants found for ', munchObj.input)
            listDF.append(tmpDF)
    return pd.concat(listDF, axis=0).reset_index(drop=True)

# Define and allelic series class
# This is the format that I'll use for are the following data analysis and visualisation
class AllelicSeries:
    def __init__(self, df): #, stan_model, stan_data, stan_fit):
        self.dfRow = df
    # Parse the data to have a wide format
        self.dfParse = df[df['Unnamed: 0'].isin(['intercept', 'slope'])] \
                            .pivot(
                                index=['gene','phenotype'], 
                                columns='Unnamed: 0', 
                                values=['Mean', 'StdDev', '2.5%', '10%', '25%', '75%', '90%', '97.5%']) \
                            .reset_index()
    # check that the data is in the right format
    # i.e. it contains the right columns, such as Mean_slope, Std deviation, et
    # This function will be used to normalise the beta slope estimated by std
    # This is to make the slope ranking more suited for prioritisation
    def parseColName(self):
        self.dfParse.columns = ['_'.join(col).strip() for col in self.dfParse.columns.values]
        return self.dfParse
    def wgtdBeta(self, betaSlopeCol = 'Mean_slope', stdDevCol = 'StdDev_slope' ):
        # This is to normalise the beta slope estimated by std
        # This is to make the slope ranking more suited for prioritisation
        tmpDF = self.parseColName()
        tmpDF['InvStdDev'] = 1 / tmpDF[stdDevCol]
        sumInvStdDev = tmpDF['InvStdDev'].sum()
        tmpDF['normMeanSlope'] = tmpDF[betaSlopeCol] * tmpDF['InvStdDev'] / sumInvStdDev
        return tmpDF
    # This function will be used to sacle between 0 and 1 and normalise the estimetes
    def sclaedNormalBeta(self, betaSlopeCol = 'Mean_slope', stdDevCol = 'StdDev_slope'):
        from sklearn.preprocessing import QuantileTransformer
        from sklearn.preprocessing import MinMaxScaler
        # This is to normalise the beta slope estimated by std
        # This is to make the slope ranking more suited for prioritisation
        tmpDF = self.wgtdBeta( betaSlopeCol, stdDevCol )
        QT = QuantileTransformer(
            output_distribution='normal', 
            random_state=0)
        MM = MinMaxScaler()
        absolute = np.abs( tmpDF['normMeanSlope'] ).values.reshape(-1, 1)
        # normalised = QT.fit_transform(absolute)
        # tmpDF['priority_score'] = MM.fit_transform( normalised )
        tmpDF['priority_score'] = MM.fit_transform( absolute )
        return tmpDF
