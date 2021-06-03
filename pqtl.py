import requests
import argparse
import logging
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

# args
parser = argparse.ArgumentParser(description='Extract robust pQTLs from epigraph db')
parser.add_argument('--gene', dest='gene', required=False, help='HGNC name')
args = parser.parse_args()

# cypher query
if args.gene is None:
    q = "MATCH (e:Exposure)<-[r:INST_EXP { tier: 'Tier1', trans_cis:'cis' }]-(i:Instruments) RETURN r;"
else:
    q = "MATCH (e:Exposure)<-[r:INST_EXP { tier: 'Tier1', trans_cis:'cis' }]-(i:Instruments) WHERE e.expID = \"{}\" RETURN r;".format(args.gene)

logging.debug("Cypher query: {}".format(q))

# get json
r = requests.get("http://api.epigraphdb.org/raw_cypher/", headers={"Accept": "application/json"}, params={'query': q, 'db': "pqtl", "api_key": "6W6ALgILczGe2ONWsFFDUxvrLTDmAPSBWbIn04Eox"})
assert r.status_code == 200

# get output
results = dict()
for r in r.json()['results']:

    if r['r']['expID'] not in results:
        results[r['r']['expID']] = []
    
    results[r['r']['expID']].append(r['r'])

# write out to file
for gene in results:

    if ';' in gene:
        logging.info("Skipping {} because the probe could not uniquely identify the protein".format(gene))
        continue

    with open(gene + ".pqtl.txt", "w") as f:
        f.write("SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\n")

        for snp in results[gene]:
            
            # get chromosome position
            variant_req = requests.get("https://grch37.rest.ensembl.org/variation/homo_sapiens/" + snp['rs_ID'], headers={ "Content-Type" : "application/json", "Accept" : "application/json"})
            assert variant_req.status_code == 200
            anno = variant_req.json()
            mapping = list(filter(lambda x: (x['assembly_name'] == "GRCh37" and x['coord_system'] == "chromosome"), anno['mappings']))[0]

            if len(mapping) == 0:
                logging.warning("Skipping {}; could not get chromosome mapping".format(snp['rs_ID']))
                continue
    
            # write to file
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(snp['rs_ID'], mapping['seq_region_name'], mapping['start'], snp['ea'].upper(), snp['nea'].upper(), snp['beta_exp'], snp['se_exp'], snp['pval_exp']))
