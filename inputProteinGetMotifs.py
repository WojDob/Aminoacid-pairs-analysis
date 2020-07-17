import json 
#AT1G65440.1  AT5G04290.1  AT2G40030.1

query = ["AT5G04290.1","AT1G65440.1","AT2G40030.1","AT3G51940.2"]

with open("allSignificantProteins.json", 'r') as f:
    significantProteins = json.load(f)

for ident in query:
	print("\n " + ident)
	for protein in significantProteins:

		if protein["id"] == ident:
			print("m {} c {} p {} z {}".format(protein["motif"],protein["count"],protein["pvalue"],protein["zscore"]))