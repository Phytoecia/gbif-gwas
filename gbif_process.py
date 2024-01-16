# %% [markdown]
# This code is to calculate the herbivore density over plant samples based on two distribution input datasets.
# It gets plant sampling data (e.g., 1001 genome accession coordinate) and GBIF data as input.
# The output file contains plant ID and corresponding density of herbivore.

# %%
import pandas as pd
import numpy as np
import geopy
from geopy import distance #distance extraction
import argparse

# %%
#default setting
query_coord='/home/minsoo/240108_GBIF_GWAS/raw_data/1001genome_coord.csv'
# target_dir='/home/minsoo/240108_GBIF_GWAS/raw_data/target_list.txt'
output_dir='./gbif_output.pheno'
intensity=True
range=100


# %%

parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", default=query_coord) #directory of arabidopsis accession and coordinate
parser.add_argument("-t", "--target", required=True) #txt file containing (gbif output dir, species name)
parser.add_argument("-o", "--outdir", default=output_dir) #output file name
parser.add_argument("-r", "--range", default=100) # range to include gbif data
parser.add_argument('-i', "--intensity", default=True)
args = parser.parse_args()


# %%
def ProcessGBIF(dir):
    raw_file = pd.read_csv(dir, sep='\t', low_memory=False)
    file_woNA = raw_file.dropna(subset=['decimalLatitude', 'decimalLongitude'])
    return file_woNA

# %%
def GetHerbivory(plant_dir, target_list, outdir, range, intensity):
    print(f"range is {range}")
    f=pd.read_csv(target_list, header=None)
    plant=pd.read_csv(plant_dir)
    plant=plant.dropna(subset=['Latitude', 'Longitude'])
    
    output=pd.DataFrame()
    output['FID']=plant['ID']
    output['IID']=plant['ID']
    
    for i, herbivore in f.iterrows():
        print(f"{herbivore[1]} calculation started")    
        density_list=[]
        for idx, ecotype in plant.iterrows():
            Herbivore_count=0
            At_coord=(ecotype['Latitude'], ecotype['Longitude'])
            for i, row in ProcessGBIF(herbivore[0]).iterrows():
                Herbivore_coord=(row['decimalLatitude'], row['decimalLongitude'])
                if geopy.distance.distance(At_coord, Herbivore_coord).km<range:
                    Herbivore_count+=1
                    if intensity==False:
                        break
            density_list.append(Herbivore_count)
            if idx%100==0:
                print(f"{idx} per 1134 accessions finished")
        output[herbivore[1]]=density_list
        print(f"{herbivore[1]} calculation finished")
    
    output.to_csv(outdir, sep=' ', index=False)
    
        

# %%
GetHerbivory(plant_dir=args.query, target_list=args.target, outdir=args.outdir, range=args.range, intensity=args.intensity)


