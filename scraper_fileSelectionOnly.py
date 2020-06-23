import sys
from urllib.parse import urlparse
import httplib2 as http
import json
import pandas as pd
from Bio import SeqIO
import urllib
import tkinter as tk
from tkinter import filedialog
import os, os.path
from tqdm import tqdm
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
import requests

'''
example with uniprot_id: cg03278611
'''

'''
Using the original file
Taking chunks of 100000 records
Using API's as headers, with json.
'''

chunksize = 1000000  # set to 1 million, just take all rows in a file (=450.000)
headers = {
    'Accept': 'application/json',
}
method = 'GET'
body = ''
h = http.Http()


# This function has been marked as debugged, tested and meets the definition of done.
#
# def promptUser():
#     """
#     Check validity of probes
#     """
#     valid_ids = []
#     user_input = [item.strip() for item in input("Enter probe id(s) and separate by comma: ").split(',')]
#     for item in user_input:
#         if not item.startswith('cg') or not len(item) == 10:
#             print(' Some IDs appears to be invalid. Check your input thoroughly. | Status: Terminated.')
#             sys.exit()
#         if item.startswith('cg') and len(item) == 10:
#             valid_ids.append(item)
#
#
#     return valid_ids

def readCGfromFiles(methylation_files):
    cgs = []
    for file in methylation_files:
        cg = pd.read_csv(file,  # not really relevant but protects against errors such as mixed types
                         low_memory=False,  # ignores DtypeWarnings
                         error_bad_lines=False,
                         skiprows=7)
        cg = cg['IlmnID'].values
        cgs.append(cg)

    ids = [item for sublist in cgs for item in sublist][:10] # The amount of records you want to run. 

    return ids



def selectedFiles():
    methylation_files = []
    print('SELECT DIRECTORY WHERE METHYLATION FILES ARE LOCATED (CHECK YOUR DESKTOP SCREEN FOR UI)')
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askdirectory()
    print('Path selected:', file_path)

    if not os.path.isdir(file_path):
        print("The directory appears to be invalid")
    if os.path.isdir(file_path):
        for file in os.listdir(file_path):
            if file.startswith('HumanMethylation') and file.endswith(".csv"):
                dir_and_file = os.path.join(file_path + '/' + file)
                methylation_files.append(dir_and_file)

    return methylation_files


# This function has been marked as debugged, tested and meets the definition of done.
def searchIDinDataframe(sub_humanMethylation_df, ids):
    """
    Get probe from dataframe
    """

    ilmn_ID = sub_humanMethylation_df[sub_humanMethylation_df['IlmnID'].isin(ids)]

    return ilmn_ID


# This function has been marked as debugged, tested and meets the definition of done.
def readDataFrame(filenames, ids):
    """
    Reads dataframe, appends probe in dataframe
    """
    appended_chunks = []
    for filename in filenames:
        for sub_humanMethylation_df in pd.read_csv(filename,
                                                   chunksize=chunksize,
                                                   # not really relevant but protects against errors such as mixed types
                                                   low_memory=False,  # ignores DtypeWarnings
                                                   error_bad_lines=False,
                                                   skiprows=7,
                                                   usecols=['IlmnID', 'Chromosome_36', 'Coordinate_36', 'Strand',
                                                            'Probe_SNPs', 'UCSC_RefGene_Name',
                                                            'UCSC_RefGene_Group', 'UCSC_CpG_Islands_Name',
                                                            'Relation_to_UCSC_CpG_Island', 'Enhancer']):
            ilmn_ID = searchIDinDataframe(sub_humanMethylation_df, ids)
            appended_chunks.append(ilmn_ID)

    df = pd.concat(appended_chunks)
    return df


# This function has been marked as debugged, tested and meets the definition of done.
def extractHugoIDs(output_probes):
    """
    Uses UCSC_RefGene_Name to lookup for genes in Hugo database
    """
    refgene_names = []
    ugrn = output_probes[['IlmnID', 'UCSC_RefGene_Name']]
    for index, row in ugrn.iterrows():
        id = row['IlmnID']
        refgene_name = row['UCSC_RefGene_Name']
        if isinstance(refgene_name, str):
            refgene_names.append([id, refgene_name])
        else:
            print(
                'extractHugoIDs() - RefGene_name not an Instance of type. (Most likely a NAN value) | status: handled')

    cg_hugo = []
    for item in refgene_names:
        cg = item[0]
        hugo_code_unique = item[1]
        if ";" in hugo_code_unique:
            items = hugo_code_unique.split(";")
            cg_hugo.append([cg, list(set(items))])
        else:
            cg_hugo.append([cg, [hugo_code_unique]])

    return cg_hugo


# This function has been marked as debugged, tested and meets the definition of done.
def getHugo(headers, method, body, h, id_uniprot_items_list):  
    """
    Gets attributes from Hugo database, with API of Hugo
    Hugo is an xml file
    return: list
    """
    hugo_record_list = []
    for record in tqdm(id_uniprot_items_list):
        id = record[0]
        for uid in record[1]:
            url = 'http://rest.genenames.org/fetch/symbol/' + uid
            target = urlparse(url)

            response, content = h.request(
                target.geturl(),
                method,
                body,
                headers)

            if response['status'] == '200':
                try:
                    data = json.loads(content)
                    name = data['response']['docs'][0]['name']
                    symbol = data['response']['docs'][0]['symbol']
                    locus_type = data['response']['docs'][0]['locus_type']
                    uniprot_ids = data['response']['docs'][0]['uniprot_ids']
                    hugo_record_list.append([id, name, symbol, locus_type, uniprot_ids])

                except KeyError as err:
                    str_err = str(err)[1:-1]  # removes quotes from error logging
                    if 'name' == str_err:
                        print('getHugo()', '- Missing value:', str_err, 'for hugo identifier:', uid,
                              '| status: handled.')
                        name = 'NA'
                        hugo_record_list.append([id, name, symbol, locus_type, uniprot_ids])
                    if 'symbol' == str_err:
                        print('getHugo()', '- Missing value:', str_err, 'for hugo identifier:', uid,
                              '| status: handled.')
                        symbol = 'NA'
                        hugo_record_list.append([id, name, symbol, locus_type, uniprot_ids])
                    if 'locus_type' == str_err:
                        print('getHugo()', '- Missing value:', str_err, 'for hugo identifier:', uid,
                              '| status: handled.')
                        locus_type = 'NA'
                        hugo_record_list.append([id, name, symbol, locus_type, uniprot_ids])
                    if 'uniprot_ids' == str_err:
                        print('getHugo()', '- Missing value:', str_err, 'for hugo identifier:', uid,
                              '| status: handled.')
                        uniprot_ids = 'NA' 
                        hugo_record_list.append([id, name, symbol, locus_type, uniprot_ids])
                    else:
                        print('ok.')

                except IndexError:
                    print('getHugo() - Hugo XML completely empty for hugo identifier:', uid, '| status: handled.')
                    continue

            else:
                print('Error detected: ' + response['status'])

    """
    used probes: cg00035864, cg00121626, cg00214611, cg00063477, cg00212031, cg01086462
    """
    return hugo_record_list


def getUniprot(hugo_data):
    """
    Parses Uniprot in xml format with SeqIO
    return: list
    """
    uniprot_records = []
    for lst in tqdm(hugo_data):
        for element in lst:
            if isinstance(element, list):
                uid = element[-1]
                if uid == 'NA':
                    print('getUniprot() - Missing value: uniprot identifier:', lst[0], '| status: handled.')
                else:
                    # XML
                    cg_id = lst[0]
                    url = 'https://www.uniprot.org/uniprot/' + uid + '.xml'
                    handle = urllib.request.urlopen(url)
                    #handle = urllib.request.urlopen(url)
                    record = SeqIO.read(handle, "uniprot-xml")
                    anno = record.annotations

                    function = anno.get('comment_function')
                    subcellular_location = anno.get('comment_subcellularlocation_location')
                    tissue_specificity = anno.get('comment_tissuespecificity')
                    disease_summary = anno.get('comment_disease')

                    if function is None:
                        function = 'NA'
                    if subcellular_location is None:
                        subcellular_location = 'NA'
                    if tissue_specificity is None:
                        tissue_specificity = 'NA'
                    if disease_summary is None:
                        disease_summary = 'NA'

                    # TXT
                    url_txt = 'http://www.uniprot.org/uniprot/' + uid + '.txt'
                    req = urllib.request.Request(url_txt)
                    page = urllib.request.urlopen(req)
                    data = page.read().decode("utf8")
                    lines = data.splitlines()

                    cell_component = ""  
                    molecular_function = ""
                    bio_process = ""
                    pathway = ""
                    
                    

                    for line in lines:
                        if " C:" in line:
                            cell_component += line[line.index(" C:") + 3:]
                        if " F:" in line:
                            molecular_function += line[line.index(" F:") + 3:]
                        if " P:" in line:
                            bio_process += line[line.index(" P:") + 3:]
                        if " Reactome;" in line:
                            pathway += line[line.index(" Reactome;") + 0:]
                            
                    cell_component = cell_component.replace('.', ' ').strip(';')
                    molecular_function = molecular_function.replace('.', ' ')
                    bio_process = bio_process.replace('.', ' ')
                    pathway = pathway.replace('.', ' ')
                    
                    if cell_component is None:
                        cell_component = 'NA'
                    if molecular_function is None:
                        molecular_function = 'NA'
                    if bio_process is None:
                        bio_process = 'NA'
                    if pathway is None:
                        pathway = 'NA'

                    uniprot_records.append([cg_id, function, subcellular_location,
                                            tissue_specificity, disease_summary,
                                            cell_component, molecular_function,
                                            bio_process, pathway, url[:-4]])

    return uniprot_records


def mergeData(output_probes, hugo_data, uniprot_data):
    """
    converts lists to dataframe:
        getHugo and getUniprot
    merges dataframes:
        hugo_df, uniprot_df and output_probes
    return: df
    """
    index = 'IlmnID'
    hugo_df = pd.DataFrame.from_records(hugo_data, columns=['IlmnID', 'name', 'symbol', 'locus_type', 'uniprot_ids'])
    uniprot_df = pd.DataFrame(uniprot_data, columns=['IlmnID', 'function', 'subcellular_location', 'tissue_specificity',
                                                     'disease_summary',
                                                     'cell_component', 'molecular_function', 'biological_process',
                                                     'pathway', 'source'])
    output_probes_df = output_probes.set_index(index)
    hugo_df = hugo_df.set_index(index)
    uniprot_df = uniprot_df.set_index(index)
    df = output_probes_df.join(hugo_df).join(uniprot_df)
    df = df.reset_index()
    df.drop_duplicates(subset='IlmnID', keep='first', inplace=True)

    return df


# ids = promptUser()

graphviz = GraphvizOutput()
graphviz.output_file = 'basic.jpg'
with PyCallGraph(output=graphviz):
    
    methylation_files = selectedFiles()
    ids = readCGfromFiles(methylation_files)
    output_probes = readDataFrame(methylation_files, ids)
    id_uniprot_items_list = extractHugoIDs(output_probes)
    hugo_data = getHugo(headers, method, body, h, id_uniprot_items_list)
    uniprot_data = getUniprot(hugo_data)
    result = mergeData(output_probes, hugo_data, uniprot_data)
    result.fillna("NA", inplace=True)
    path_out = 'C:\\Users\\aditi\\OneDrive\\Bureaublad\\webscraper- adi\\main\\webvisualisation\\data'
    result.to_csv(path_out, "result.csv")
    print('Done.')
