import pandas as pd
import sys
from munch import Munch

def parse_vep(vep_file):
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

            try:
                transcriptConsequencesDF = pd.json_normalize(munchObj.transcript_consequences) # This is all parsed correctly in this way
            except(KeyError, AttributeError):
                # go to next line if no transcript consequences
                print( munchObj.input, 'doesn t have canonical transcript')
                continue
            if transcriptConsequencesDF.columns.str.contains('canonical').any():
                pass
            else:
                print( munchObj.input, 'doesn t have canonical transcript')
                continue
            
            transcriptConsequencesDF = transcriptConsequencesDF[transcriptConsequencesDF.canonical == 1]
            
            # get the line with hte list missing information
            transcriptConsequencesDF['missingInfo'] = transcriptConsequencesDF.apply(lambda x: x.isna().sum(), axis=1)
            transcriptConsequencesDF = transcriptConsequencesDF[ transcriptConsequencesDF.missingInfo == transcriptConsequencesDF.missingInfo.min()]
            
            colKeep = ['gene_id','primateai', 'revel', 'cadd_phred', 'cadd_raw', 'pli_gene_value', 
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

read_files = []
for file in sys.argv[1:]:
    parsedFile = parse_vep(file)
    read_files.append(parsedFile)

pd.concat(read_files, axis=0) \
    .reset_index(drop=True) \
    .to_csv('parsed_vep.csv', 
            index=False)