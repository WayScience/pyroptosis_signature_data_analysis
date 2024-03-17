#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pathlib

import pandas as pd

# In[3]:


# set the ouput path for the library file
output_path = pathlib.Path("../IDR_metadata/screenA/idr0000-screenA-library.txt")


# In[5]:


# set path to the metadata file
path = pathlib.Path("../../data/Interstellar_plate2_platemap.csv")

# read the metadata file
df = pd.read_csv(path)

# remove any rows that have NaN in the cell_type column
df = df.dropna(subset=["cell_type"])

# add a Plate column to the dataframe
df["Plate"] = "plate_1"
# rename the columns
df = df.rename(columns={"well_id": "Well", "cell_type": "Characteristics[cell type]"})
df["Term Source 1 REF"] = "NCBITaxon"
df["Term Source 1 Accession"] = "NCBITaxon_9606"
df.drop(columns="incubation inducer (h)", inplace=True)
df[
    "Channels"
] = "Hoechst 33342:nucleus;concanavalin A (con A) PhenoFluor488 conjugate:endoplasmic reticulumn;wheat germ agglutinin (WGA) PhenoFluor555 conjugate:Golgi apparatus and plasma membrane;phalloidin PhenoFluor568 conjugate: F-actin;Mitochondrial stain PhenoVue 641: mitochondria; AlexaFluor514: Cleaved N-terminal Gasdermin D"
df.head()


# In[6]:


# create a dictionary to map the compound names to the compound IDs
compund_dicts = {
    "IUPAC_names": {
        "DMSO": "methylsulfinylmethane",
        "Disulfiram": "diethylcarbamothioylsulfanyl N,N-diethylcarbamodithioate",
        "Topotecan": "(19S)-8-[(dimethylamino)methyl]-19-ethyl-7,19-dihydroxy-17-oxa-3,13-diazapentacyclo[11.8.0.02,11.04,9.015,20]henicosa-1(21),2,4(9),5,7,10,15(20)-heptaene-14,18-dione",
        "H202": "hydrogen peroxide",
        "Thapsigargin": "[(3S,3aR,4S,6S,6aR,7S,8S,9bS)-6-acetyloxy-4-butanoyloxy-3,3a-dihydroxy-3,6,9-trimethyl-8-[(Z)-2-methylbut-2-enoyl]oxy-2-oxo-4,5,6a,7,8,9b-hexahydroazuleno[4,5-b]furan-7-yl] octanoate",
        "Nigericin": "sodium;(2R)-2-[(2R,3S,6R)-6-[[(2S,4R,5R,6R,7R,9R)-2-[(2R,5S)-5-[(2R,3S,5R)-5-[(2S,3S,5R,6R)-6-hydroxy-6-(hydroxymethyl)-3,5-dimethyloxan-2-yl]-3-methyloxolan-2-yl]-5-methyloxolan-2-yl]-7-methoxy-2,4,6-trimethyl-1,10-dioxaspiro[4.5]decan-9-yl]methyl]-3-methyloxan-2-yl]propanoate",
        "Z-VAD-FMK": "methyl (3S)-5-fluoro-3-[[(2S)-2-[[(2S)-3-methyl-2-(phenylmethoxycarbonylamino)butanoyl]amino]propanoyl]amino]-4-oxopentanoate",
    },
    "InCHI_keys": {
        "DMSO": "IAZDPXIOMUYVGZ-UHFFFAOYSA-N",
        "Disulfiram": "AUZONCFQVSMFAP-UHFFFAOYSA-N",
        "Topotecan": "UCFGDBYHRUNTLO-QHCPKHFHSA-N",
        "H202": "MHAJPDPJQMAIIY-UHFFFAOYSA-N",
        "Thapsigargin": "IXFPJGBNCFXKPI-FSIHEZPISA-N",
        "Nigericin": "MOYOTUKECQMGHE-PDEFJWSRSA-M",
        "Z-VAD-FMK": "MIFGOLAMNLSLGH-QOKNQOGYSA-N",
    },
    "SMILES": {
        "DMSO": "CS(=O)C",
        "Disulfiram": "CCN(CC)C(=S)SSC(=S)N(CC)CC",
        "Topotecan": "CCC1(C2=C(COC1=O)C(=O)N3CC4=CC5=C(C=CC(=C5CN(C)C)O)N=C4C3=C2)O",
        "H202": "OO",
        "Thapsigargin": "CCCCCCCC(=O)OC1C2C(=C(C1OC(=O)C(=CC)C)C)C3C(C(CC2(C)OC(=O)C)OC(=O)CCC)(C(C(=O)O3)(C)O)O",
        "Nigericin": "CC1CCC(OC1C(C)C(=O)[O-])CC2CC(C(C3(O2)C(CC(O3)(C)C4CCC(O4)(C)C5C(CC(O5)C6C(CC(C(O6)(CO)O)C)C)C)C)C)OC.[Na+]",
        "Z-VAD-FMK": "CC(C)C(C(=O)NC(C)C(=O)NC(CC(=O)OC)C(=O)CF)NC(=O)OCC1=CC=CC=C1",
    },
    "PubChem_urls": {
        "DMSO": "https://pubchem.ncbi.nlm.nih.gov/compound/679/",
        "Disulfiram": "https://pubchem.ncbi.nlm.nih.gov/compound/3117/",
        "Topotecan": "https://pubchem.ncbi.nlm.nih.gov/compound/60700/",
        "H202": "https://pubchem.ncbi.nlm.nih.gov/compound/784/",
        "Thapsigargin": "https://pubchem.ncbi.nlm.nih.gov/compound/446378",
        "Nigercin": "https://pubchem.ncbi.nlm.nih.gov/compound/16760591",
        "Z-VAD-FMK": "https://pubchem.ncbi.nlm.nih.gov/compound/5497174/",
    },
    "Control_type": {
        "DMSO": "Negative",
        "Disulfiram": "",
        "Topotecan": "",
        "H202": "",
        "Thapsigargin": "",
        "Nigericin": "",
        "Z-VAD-FMK": "",
    },
}


# In[7]:


# add a column with values mapped from the compound_dict
df["Compound1"] = df["inducer1"].map(compund_dicts["IUPAC_names"])
df["Compound1_InChIKey"] = df["inducer1"].map(compund_dicts["InCHI_keys"])
df["Compound1_SMILES"] = df["inducer1"].map(compund_dicts["SMILES"])
df["Compound1_PubChem_URL"] = df["inducer1"].map(compund_dicts["PubChem_urls"])
df["Compound1_Control_type"] = df["inducer1"].map(compund_dicts["Control_type"])

df["Compound2"] = df["inducer2"].map(compund_dicts["IUPAC_names"])
df["Compound2_InChIKey"] = df["inducer2"].map(compund_dicts["InCHI_keys"])
df["Compound2_SMILES"] = df["inducer2"].map(compund_dicts["SMILES"])
df["Compound2_PubChem_URL"] = df["inducer2"].map(compund_dicts["PubChem_urls"])
df["Compound2_Control_type"] = df["inducer2"].map(compund_dicts["Control_type"])

df["Compound3"] = df["inhibitor"].map(compund_dicts["IUPAC_names"])
df["Compound3_InChIKey"] = df["inhibitor"].map(compund_dicts["InCHI_keys"])
df["Compound3_SMILES"] = df["inhibitor"].map(compund_dicts["SMILES"])
df["Compound3_PubChem_URL"] = df["inhibitor"].map(compund_dicts["PubChem_urls"])
df["Compound3_Control_type"] = df["inhibitor"].map(compund_dicts["Control_type"])


# In[8]:


# replace all NaN values with an empty string
df = df.fillna("")
df.head()


# In[9]:


# save the dataframe to a tab-delimited file
df.to_csv(output_path, sep="\t", index=False)
