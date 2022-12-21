import argparse
import pantherAutomation as pa
from os.path import split, splitext, abspath, exists, join
from os import mkdir
from datetime import date



def main():
	args, exp, back = read_arguments()
	GO_response = run_panther_overexpression_test(args, exp, back)
	GO_signif, all_GO = pa.filter_GO_results(GO_response)
	pa.save(GO_signif, outpath = args.outputdir, experimentName=args.experimentName, annot=args.annotation, signif=True)
	pa.save(all_GO, outpath = args.outputdir, experimentName=args.experimentName, annot=args.annotation, signif=False)


def read_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', type=str, 
		required=True, help='Path to input file containing genes of interest')
	parser.add_argument('-r', '--reference', type=str, 
		required=True, help='Path to reference (background) gene list')
	parser.add_argument('-o', '--organism', type=int, default=10090, 
		required=False, help='A numeric ID representing the organism')
	parser.add_argument('-a', '--annotation', type=str, default='GO:0008150', 
		required=False, help='The ontology(ies) in which you want to look for enrichment.\
		Default is GO-BP. Additional options are GO:0003674 (GO-MF), and GO:0005575 (GO-CC)')
	parser.add_argument('-s', '--stats', type=str, default='FISHER',required=False, 
		choices=['FISHER', 'BINOMIAL'], help='The statistical test performed on the enrichment values')
	parser.add_argument('-e', '--experimentName', type=str, default=f'Exp_{date.today()}', required=False)
	parser.add_argument('-c', '-correction', type=str, default='FDR', required=False,
		choices=['FDR', 'BONFERRONI', 'NONE'], help='Multiple testing pvalue correction method')
	parser.add_argument('-z', '--outputdir', required=True, type=str)
	args = parser.parse_args()
	print(args.annotation)
	# Generate the comma separated background and comparison lists 
	backgroundList = pa.parse_gene_list(args.reference)
	experimentList = pa.parse_gene_list(args.input)

	# In outdir, make two folders, if they aren't already there for all results and significant results
	myoutpath = abspath(args.outputdir)
	new_folders = ['significant_pantherGO_results', 'all_pantherGO_results']
	dirs = [join(myoutpath, fldr) for fldr in new_folders]
	for newdir in dirs:
		if not exists(newdir):
			mkdir(newdir)
		else:
			continue
	return args, backgroundList, experimentList

def run_panther_overexpression_test(args, experimentList, backgroundList):
	parameters: dict = pa.package_request_parameters(inputList=experimentList, backgroundList=backgroundList, 
		refListorganism=args.organism, annotDataSet=args.annotation, enrichmentTestType=args.stats, 
		correction='FDR')
	base_url = 'http://pantherdb.org/services/oai/pantherdb/enrich/overrep?'
	query_response_json = pa.call_request(base_url, parameters)
	filename: str = splitext(split(args.input)[1])[0]
	response_df = pa.reformat_response_to_df(query_response_json)
	return response_df


if __name__ == '__main__':
	main()




