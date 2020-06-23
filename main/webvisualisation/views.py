import inspect
import os
import pandas as pd
from django.shortcuts import render
from .forms import InputForm
import csv


download_folder = os.path.join('static') + '/downloads/'
print(download_folder)

def export(vars, url):
    with open(url, 'w', newline='') as csvfile:
        record_writer = csv.writer(csvfile, delimiter=',')
        record_writer.writerow(vars)

def index(request):
    form = InputForm()
    if request.method == 'GET':
        context = {'form': form}
        return render(request, 'index.html', context)

    if request.method == 'POST':
        form = InputForm(request.POST)
        if form.is_valid():
            cg_code = form.cleaned_data['cg_code']
            directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/data/result.csv'
            directory = directory.replace('\\', '/')
            df = pd.read_csv(directory)
            df.reset_index(inplace=True)
            row = df.loc[df['IlmnID'] == cg_code]

            if row.empty:
                IlmnID = 'No Data'
                name = 'No Data'
                symbol = 'No Data'
                locus_type = 'No Data'
                uniprot_ids = 'No Data'
                function = 'No Data'
                subcellular_location = 'No Data'
                tissue_specificity = 'No Data'
                disease_summary = 'No Data'
                cell_component = 'No Data'
                molecular_function = 'No Data'
                biological_process = 'No Data'
                pathway = 'No Data'
                source = 'No Data'

            else:
                IlmnID = row['IlmnID'].values.item()
                name = row['name'].values.item()
                symbol = row['symbol'].values.item()
                locus_type = row['locus_type'].values.item()
                uniprot_ids = row['uniprot_ids'].values.item()
                function = row['function'].values.item()
                subcellular_location = row['subcellular_location'].values.item()
                tissue_specificity = row['tissue_specificity'].values.item()
                disease_summary = row['disease_summary'].values.item()
                cell_component = row['cell_component'].values.item()
                molecular_function = row['molecular_function'].values.item()
                biological_process = row['biological_process'].values.item()
                pathway = row['pathway'].values.item()
                source = row['source'].values.item()

            vars = [IlmnID, name, symbol, locus_type, uniprot_ids, function, subcellular_location,
                    tissue_specificity, disease_summary, cell_component, molecular_function,
                    biological_process, pathway, source]

            url = download_folder + IlmnID + '.csv'

            export(vars, url)

            context = {'form': form, 'IlmnID': IlmnID, 'url': url, 'name': name, 'symbol': symbol, 'locus_type': locus_type,
                       'uniprot_ids': uniprot_ids, 'function': function,
                       'subcellular_location': subcellular_location,
                       'tissue_specificity': tissue_specificity, 'disease_summary': disease_summary,
                       'cell_component': cell_component, 'molecular_function': molecular_function,
                       'biological_process': biological_process, 'pathway': pathway, 'source': source}

            return render(request, 'index.html', context)

        else:
            context = {'form': form}
            return render(request, 'index.html', context)

    context = {'form': form}
    return render(request, 'index.html', context)
