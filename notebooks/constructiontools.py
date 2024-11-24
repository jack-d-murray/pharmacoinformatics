# Note: add sleep to ensure compliance with current request limits where appropriate. 
from chembl_structure_pipeline import standardizer
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.MolStandardize import rdMolStandardize

# Data manipulation libraries
import numpy as np
import pandas as pd

# Web scraping libraries
from bs4 import BeautifulSoup
from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings
import pubchempy as pcp
import requests

# Other necessary libraries
import ast
import time
from collections import Counter
from functools import reduce
from math import gcd
import os
from pdfminer.high_level import extract_text
import regex as re

def add_excipients(register):
    """Updates the split register to add the excipients raw regular expression match from the split txt file
    
    parameters: EMA Central Register split into products (pandas DataFrame)
    return: The above dataframe with the extra column, excipients_raw (pandas DataFrame)
    """
    register["excipients_raw"] = " "
    pattern = re.compile(r"6.1.?\s*List of excipients((.|\n)*)6\.2.?\s*Incom")
    for index, row in register.iterrows():
        try:
            document = register.at[index, "Document"]
            with open(f"./regulatory_docs/txt_docs_split/{document}.txt", encoding="utf-8") as f:
                text = f.read()
                matches = pattern.findall(text)
                list_of_excipients = matches[0][0].strip().upper().replace("\n", ",").split(",")
                list_of_excipients_stripped = [excipient.strip() for excipient in list_of_excipients if excipient.strip()]
                register.at[index, "excipients_raw"] = list_of_excipients_stripped
        except:
            register.at[index, "excipients_raw"] = "Unsuccessful"
    return register


def add_pharmaceutical_form(register):
    """Updates the split register to add the pharmaceutical form raw regular expression match from the split txt file
    
    parameters: EMA Central Register split into products (pandas DataFrame)
    return: The above dataframe with the extra column, pharmaceutical_form (pandas DataFrame)
    """
    register["pharmaceutical_form"] = " "
    pattern = re.compile(r"(?i)PHARMACEUTICAL\s?\s?\s?FORM((.|\n)*)4.?\s*CLINI?CAL\s?\s?\s?PARTICULARS")
    for index, row in register.iterrows():
        try:
            document = register.at[index, "Document"]
            with open(f"../regulatory_docs/txt_docs_split/{document}.txt", encoding="utf-8") as f:
                text = f.read()
                matches = pattern.findall(text)
                register.at[index, "pharmaceutical_form"] = matches[0][0][:10000].lower().strip().replace(",", "-")
        except:
            register.at[index, "pharmaceutical_form"] = "Unsuccessful"
    return register

def add_qq_composition(register):
    """Updates the split register to add the quantitative and qualitative composition from the split txt file
    
    parameters: EMA Central Register split into products (pandas DataFrame)
    return: The above dataframe with the extra column, qq_composition (pandas DataFrame)
    """
    register["qq_composition"] = " "
    for index, row in register.iterrows():
        document = register.at[index, "Document"]
        with open(f"../regulatory_docs/txt_docs_split/{document}.txt", encoding="utf-8") as f:
            text = f.read()
            register.at[index, "qq_composition"] = text[:400].lower().strip().replace(",", "-")
    return register

def add_pharmacokinetics(register):
    """Updates the split register to add the pharmacokinetics section from the split txt file
    
    parameters: EMA Central Register split into products (pandas DataFrame)
    return: The above dataframe with the extra column, pharmacokinetics (pandas DataFrame)
    """
    register["pharmacokinetics"] = " "
    pattern = re.compile(r"(?i)PHARMACOKINETIC\s?\s?\s?PROPERTIES((.|\n)*)5.?.?.?\s*PRE-?CLINICAL\s?\s?\s?SAFETY")
    for index, row in register.iterrows():
        try:
            document = register.at[index, "Document"]
            with open(f"./regulatory_docs/txt_docs_split/{document}.txt", encoding="utf-8") as f:
                text = f.read()
                matches = pattern.findall(text)
                register.at[index, "pharmacokinetics"] = matches[0][0][:10000].lower().strip().replace(",", "-")
        except:
            register.at[index, "pharmacokinetics"] = "Unsuccessful"
    return register

def remove_products(df):
    """Removes unnessesary products from raw output
    
    parameters: EMA Central Register (pandas DataFrame)
    return: EMA Central Register without veterinary, biosimilars, vaccines, and refused products (pandas DataFrame)
    """
    df = df[df["Category"] == "Human"]
    df = df[df["Authorisation status"] != "Refused"]
    df = df[df["Biosimilar"] == "no"]
    df = df[df["Human pharmacotherapeutic group"] != "Vaccines"]
    df.drop(columns = ["Biosimilar"], inplace=True)
    df.dropna(axis=1, how="all", inplace=True)
    return df

def identify_small_molecules(register, sms):
    """Removes products that do not contain small molecule drugs
    
    parameters: EMA Central Register (pandas DataFrame), Substances Management Services export (pandas Dataframe)
    return: EMA Central Register containing only small molecules (pandas Dataframe)
    """
    sms = sms[sms["Language"] == "English"]
    sms = sms[sms["Substance_Domain"] == "Human use"]
    sms["Substance_Name"] = sms["Substance_Name"].str.upper()
    sms = sms[~sms.Substance_Type.isin(["Chemical", "Specified Substance Group 1", "Specified Substance Group 3", "Chemical", "Mixture", "Nucleic acid"])]
    large_molecules_list = sms["Substance_Name"].to_list()
    register["Active substance"].fillna("Blank", inplace=True)
    register["Active substance"] = register["Active substance"].str.upper().str.split(",").apply(lambda x: [item.strip() for item in x])
    register["Large_mol"] = register["Active substance"].apply(lambda x: any(drug in x for drug in large_molecules_list))
    register.drop(index=register[register["Large_mol"]].index, inplace=True)
    register.pop("Large_mol")
    return register

def identify_nucleic_acids(register, sms):
    """Identifies molecules which are nucleic acids to make manual removal of large nucleic acids easier
    
    All nucleic acids identified by this function were manually verified to be large molecules
    
    parameters: EMA Central Register (pandas DataFrame), Substances Management Services export (pandas DataFrame)
    return: List of nucleic acids (python list)
    """
    sms = sms[sms["Language"] == "English"]
    sms = sms[sms["Substance_Domain"] == "Human use"]
    sms["Substance_Name"] = sms["Substance_Name"].str.upper()
    sms = sms[sms.Substance_Type.isin(["Nucleic acid"])]
    nucleic_acid_sms = sms["Substance_Name"].to_list()
    register["Nucleic_acid"] = register["Active substance"].apply(lambda x: any(drug in x for drug in nucleic_acid_sms))
    register.drop(index=register[register["Nucleic_acid"]].index, inplace=True)
    register.pop("Nucleic_acid")
    return register

def assign_inchi_key(register, sms):
    """Assigns an InChI Key from the Substances Management Services for each product in the central register
    
    parameters: EMA Central Register (pandas DataFrame), Substances Management Services export (pandas DataFrame)
    return: EMA Central Register with small molecules mapped to an InChI Key (pandas DataFrame)
    """
    sms = sms[sms["Language"] == "English"]
    sms = sms[sms["Substance_Domain"] == "Human use"]
    sms["Substance_Name"] = sms["Substance_Name"].str.upper()
    register["InChI_Keys"] = [list() for _ in range(len(register))]
    def substance_to_inchikey(substances):
        matching_inchikeys = []
        for substance in substances:
            matching_inchikey = sms.loc[sms["Substance_Name"] == substance, "Inchikey"].tolist()
            if matching_inchikey:
                matching_inchikeys.extend(matching_inchikey)
            else:
                matching_inchikeys.append("not found")
        return matching_inchikeys
    register["InChI_Keys"] = register["Active substance"].apply(substance_to_inchikey)
    register["InChI_Keys"] = register["InChI_Keys"].apply(lambda x: ["not found"] if all(pd.isna(inchikey) for inchikey in x) else x)
    return register

def get_inchi_from_inchikey(inchi_key):
    """SMS alternative to PubChemPy, uses RESTful URLS and SMS to obtain InChI from InChI Key
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchi_key}/property/InChI/JSON"
    response = requests.get(url)
    
    if response.status_code != 200:
        return None
    inchi_value = response.json()["PropertyTable"]["Properties"][0]["InChI"]
    return inchi_value

def clean_excipients(register, mapping_df):
    """Adds a new column, "excipients_matched", to the main table.
    
    parameters: EMA Central Register containing column excipients_raw (pandas DataFrame),
                DataFrame with excipients matched to the controlled vocabulary (pandas DataFrame)
    return: The above dataframe with the extra column, excipients_matched (pandas DataFrame)
    """
    def get_cleaned_list(excipient_string):
        try:
            excipient_list = ast.literal_eval(excipient_string)
            return [mapping_df.loc[mapping_df["Value_z"] == excipient, "Value_x"].iloc[0] 
                    if excipient in mapping_df["Value_z"].values 
                    else "BLANK" 
                    for excipient in excipient_list]
        except (ValueError, SyntaxError, KeyError, IndexError):
            return "Unsuccessful"
    register["excipients_matched"] = register["excipients_raw"].apply(get_cleaned_list)
    return register
    
def get_pubchem_for_excipient(excipient_name):
    try:
        compound = pcp.get_compounds(excipient_name, "name")[0]
        print(f"{excipient_name} - success")
        time.sleep(0.2)
        return pd.Series(compound.inchi)
    except IndexError:
        print(f"{excipient_name} - no match")
        time.sleep(0.2)
        return pd.Series(None)

def assign_identifier_to_actives(actives_df, sms):
    """Assigns an identifier from the Substances Management Services for each active in actives_df (column= active_PSS)
    
    parameters: actives_df (pandas DataFrame), Substances Management Services export (pandas DataFrame)
    return: actives_df with each excipient mapped to an FDA preferred substance name
    """
    sms = sms[sms["Language"] == "English"]
    sms = sms[sms["Substance_Domain"] == "Human use"]
    sms["Substance_Name"] = sms["Substance_Name"].str.upper()
    temp_df = pd.merge(actives_df, 
                    sms[
                        ["Substance_Name","Substance_Type", "Molecular_Formula", 
                         "Molecular_Weight", "Inchikey", "Comment", 
                         "Created_Date", "Last_Updated_Date", "External_Code_XEVMPD", 
                         "External_Code_SVG", "External_Code_UNII"]],
                    how="left",
                    left_on = "active_PSS",
                    right_on = "Substance_Name")
    def alternative_matching_strategy(drug_name):
        cid_list = pcp.get_cids(f"{drug_name}")
        if cid_list != []:
            compound = pcp.Compound.from_cid(cid_list[0])
            return compound.inchikey
        else:
            return None
    temp_df.loc[temp_df["Inchikey"].isna(), "Inchikey"] = temp_df.loc[temp_df["Inchikey"].isna(), "active_PSS"].apply(alternative_matching_strategy)
    return temp_df

def download_pdfs(register, output_dir = "../regulatory_docs/pdf_docs"):
    """ Downloads .pdf files from a pdf URL and writes to an output directory, with a default provided
    This function will not be provided as part of the published version of this paper to prevent requests being sent to the EMA
    website. Sleep can be added to reduce number of requests.
    """
    for index, pdf_url in register.iterrows():
        pdf = pdf_url["PDF_URL"]
        with requests.get(pdf) as response:
            with open(f"{output_dir}/{index}.pdf", "wb") as file:
                file.write(response.content)
        print(f"{index} - downloaded successfully")
        
def get_pdf_link(site):
    """ Obtains pdf link for product given the link to the product's main page on the EMA website. Functional as of 2/1/2024.
    This function will not be provided as part of the published version of this paper to prevent requests being sent to the EMA
    website. Sleep can be added to reduce number of requests.
    """
    try:
        result = requests.get(site)
        output = BeautifulSoup(result.text, "html.parser")
        product_section = output.find_all(class_="ema-section product-info")
        english_link_regex = re.compile(r"documents/product-information/[^/]+-epar-product-information_en\.pdf")
        correct_section = product_section[0].find("a", href=english_link_regex)
        href = correct_section["href"]
        return f"https://www.ema.europa.eu{href}"
        print("Success")
    except:
        return "Source pdf URL manually"

def update_master_table(register):
    """Updates the register to account for SmPCs containing information on more than one product
    
    parameters: EMA Central Register (pandas DataFrame)
    return: EMA Central Register split based on the products in "./regulatory_docs/txt_docs_split/" (pandas DataFrame)
    """
    doclist = []
    for file in os.listdir("../regulatory_docs/txt_docs_split/"):
        file_name = file.removesuffix(".txt")
        doclist.append(file_name)
    docframe = pd.DataFrame(doclist)
    docframe.columns = ["Document"]
    docframe["Product"] = " "
    docframe["Product"] = docframe["Document"].apply(lambda i: int(i.split("_")[0]))
    return_frame = docframe.merge(register, left_on="Product", right_on="Column1")
    return_frame.sort_values("Product")
    return return_frame

def get_inchi_from_chembl(inchikey):
    """ Get Inchi from inchikey using ChEMBL,
    Ensure compliance with request limits.
    """
    try:
        molecule = new_client.molecule.get(f"{inchikey}")
        if molecule:
            return molecule["molecule_structures"]["standard_inchi"]
        return None
    except:
        return "error"

def get_isomeric_smiles_from_inchikey(inchi_key):
    """SMS alternative to PubChemPy, uses RESTful URLS and SMS to obtain Isomeric SMILES from InChI Key
    Ensure compliance with current request limits
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchi_key}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    
    if response.status_code != 200:
        return None
    isomeric_smiles = response.json()["PropertyTable"]["Properties"][0]["IsomericSMILES"]
    return isomeric_smiles

def desolvate(isomeric_smiles):
    """Generates PS (Parent-Salt)/Desolvated form form from Isomeric Smiles
    """
    
    # Define some variables
    mol = Chem.MolFromSmiles(isomeric_smiles)
    parent_molblock, _ = standardizer.get_parent_molblock(Chem.MolToMolBlock(mol))
    water_smarts = Chem.MolFromSmarts("[OH2]") # Hydrate
    ethanolate_smarts = Chem.MolFromSmarts("[CCO]")
    water_matches = mol.GetSubstructMatches(water_smarts)
    ethanolate_matches = mol.GetSubstructMatches(ethanolate_smarts)
    remover = SaltRemover(defnData = "[OH2,CCO]")

    # To catch entries which are organometallics/nonsalts/radioactive etc.
    if Chem.MolToInchiKey(Chem.MolFromMolBlock(parent_molblock)) == Chem.MolToInchiKey(mol):
        return Chem.MolToSmiles(mol), 1
    
    # For organic drug molecules which are salts but do not contain unbonded water
    if not (water_matches or ethanolate_matches):
        component_list = isomeric_smiles.split(".")
        counts = Counter(component_list)
        return Chem.MolToSmiles(mol), counts[max(counts, key=len)] # may still have multiple drug molecules
    
    # For molecules with water nonbonded to main structure
    elif water_matches or ethanolate_matches:
        stripped = remover.StripMol(mol) # Remove water with RDKit
        no_water_smiles = Chem.MolToSmiles(stripped, isomericSmiles = True) # Generate isomeric smiles
        component_list = no_water_smiles.split(".") # Split on the non-bonded character
        counts = Counter(component_list) # Generate Counter object of counts
        hcd = gcd(*counts.values()) # Calculate highest common demoninator of counts
        for key in counts:
            counts[key] //= hcd # Divide by highest common demoninator of counts
        elements = []
        for key, count in counts.items():
            elements.extend([key] * count) # Generate new smiles string
        new_smiles = ".".join(elements)
        final_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(new_smiles)) # canonical isomeric smiles
        return final_smiles, counts[max(counts, key=len)] # Counts refers to the count of the longest component
    
    # Return value for manual curation
    else:
        return "manual", "manual"
    
def get_p_form(isomeric_smiles):
    """  Get parent form of a drug substance using ChEMBL curation protocols
    """
    try:
        mol = Chem.MolFromSmiles(isomeric_smiles)
        parent_molblock, _ = standardizer.get_parent_molblock(Chem.MolToMolBlock(mol))
        parent_form = Chem.MolFromMolBlock(parent_molblock)
        return Chem.MolToSmiles(parent_form), Chem.MolToInchi(parent_form), Chem.MolToInchiKey(parent_form) # isomericSmiles and canonical both default to true
    except:
        return "error", "error", "error" 

def get_chembl_id(inchikey):
    """ Retrieves ChEMBL ID given an InChIKey.
    Ensure compliance with request limits
    """
    try:
        molecule = new_client.molecule.get(f"{inchikey}")
        if molecule:
            print("Success")
            return molecule["molecule_chembl_id"]
    except:
        print(f"{inchikey} Unsuccessful")
        return "error"

def get_inchikey_from_inchi(inchi):
    try:
        return Chem.MolToInchiKey(Chem.MolFromInchi(inchi))
    except:
        return np.nan
        
def turn_actives_into_list(string):
    return [item.strip().upper() for item in string.split(",")]

def nan_to_list(string):
    try:
        return ast.literal_eval(string.replace("nan", "None"))
    except:
        return "Error"

def get_pubchem_cid(inchikey):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    url = f"{base_url}/inchikey/{inchikey}/cids/TXT"
    
    try:
        response = requests.get(url)
        response.raise_for_status()  
        cid = response.text.strip()
        print(f"{inchikey}: success")
        return int(cid)
    except (requests.HTTPError, ValueError):
        print("error")
        return np.nan

def create_combined_data_dictionary(dataframes_with_names):
    """ Accepts list of tuples (dataframe, dataframe_name). 
    Returns: data dictionary as Pandas Dataframe.
    """
    data_dict_list = []
    
    for df, table_name in dataframes_with_names:
        data_dict = pd.DataFrame({
            'table_name': table_name,
            'Column Name': df.columns,
            'Data Type': df.dtypes.astype(str),
            'Non-Null Count': df.notnull().sum(),
            'Unique Values': df.nunique(),
            'Sample Values': df.apply(lambda x: ', '.join(x.dropna().astype(str).unique()[:3]), axis=0)
        })
        data_dict_list.append(data_dict)
    combined_data_dict = pd.concat(data_dict_list, ignore_index=True)
    return combined_data_dict
    