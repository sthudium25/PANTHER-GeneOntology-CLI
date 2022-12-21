import requests
from os.path import abspath, isfile, split
from os import devnull
import mimetypes
from urllib.parse import quote
from collections import ChainMap
import pandas as pd
from tqdm import tqdm


ontologies = {'GO:0003674':'GO-MF', 'GO:0008150':'GO-BP', 'GO:0005575':'GO-CC'}

def parse_gene_list(path: str) -> str:
	filepath = abspath(path)
	sep = ','
	if isfile(filepath):
		ftype, _ = mimetypes.guess_type(filepath)
		if ftype == 'text/plain':
			sep = '\t'
		elif ftype != 'text/csv':
			raise TypeError('background list type must be .txt or .csv')
	with open(filepath, 'r') as f:
		ref_list = [line.rstrip() for line in f.readlines()]
		gene_list_as_str = ','.join(ref_list)
	return gene_list_as_str

def package_request_parameters(inputList, backgroundList, refListorganism: int, 
	annotDataSet: str, enrichmentTestType: str, correction: str) -> dict:
	parameters = {}
	parameters['geneInputList'] = inputList
	supported_organisms = requests.get("http://pantherdb.org/services/oai/pantherdb/supportedgenomes")
	orgs_list = supported_organisms.json()['search']['output']['genomes']['genome']
	found = False
	for i in range(len(orgs_list)):
		if orgs_list[i]['taxon_id'] == refListorganism:
			found = True
	if found: 
		parameters['organism'] = str(refListorganism)
		parameters['refInputList'] = backgroundList
		parameters['refOrganism'] = str(refListorganism)
	else:
		raise AttributeError('Specified organism is not supported by Panther')

	if annotDataSet not in ontologies.keys():
		raise  AttributeError('Specified data set not available')
	else:
		parameters['annotDataSet'] = annotDataSet

	if enrichmentTestType not in ['FISHER', 'BINOMIAL']:
		raise  AttributeError('Specified statistical test not supported')
	else:
		parameters['enrichmentTestType'] = enrichmentTestType

	if correction not in ['FDR', 'BONFERRONI', 'NONE']:
		raise AttributeError('Specified multiple testing method not supported')
	else:
		parameters['correction'] = correction
	return parameters

def call_request(request_url, parameters):
	query_res = requests.post(request_url, data=parameters, stream=True)
	# tot_size_in_bytes = int(query_res.headers.get('content-length', 0))
	# block_size = 1024
	# progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
	# for data in response.iter_content(block_size):
	# 	progress_bar.update(len(data))
		
 #    progress_bar.close()
	if query_res.status_code == 200:
		print(f'{query_res.status_code}: response returned successfully.')
		return query_res.json()
	else:
		print("Bad query. The response code returned by the server was: " + str(query_res.status_code))


def reformat_response_to_df(query_json: requests.Response) -> pd.DataFrame:
	res_json = query_json['results']['result']
	for elt in res_json:
		tmp = ChainMap(elt['term'])
		del elt['term']
		for k, v in tmp.maps[0].items():
			elt[k] = v
	query_df = pd.DataFrame(res_json)
	query_df = query_df.reindex(columns=['id','label','number_in_reference', 'number_in_list', 
    	'expected', 'fold_enrichment', 'plus_minus', 'pValue', 'fdr'])
	return query_df

def filter_GO_results(response_df: pd.DataFrame):
	response_df.sort_values(by='fdr', ascending=True, inplace=True)
	signif = response_df[response_df['fdr'] < 0.05]
	all_results = response_df
	return signif, all_results
	
	
def save(output_df: pd.DataFrame, outpath: str, experimentName: str, annot: str, signif: bool):
	new_annot = ontologies[annot]
	if signif:
		output_df.to_csv(path_or_buf=f"{outpath}/significant_pantherGO_results/{experimentName}_{new_annot}_pantherGO.csv")
	else:
		output_df.to_csv(path_or_buf=f"{outpath}/all_pantherGO_results/{experimentName}_{new_annot}_pantherGO.csv")


