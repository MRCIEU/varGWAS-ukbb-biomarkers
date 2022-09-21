#!/usr/bin/env python3

# Import relevant libraries for HTTP request and JSON formatting
import requests
import json
import sys
import pandas as pd

def query_OpenTargetGenetics_variantForRSID(variant_id):
    # Build query string
    query_string = """
    query genesForVariant ($myVariantId: String!){
      search(queryString: $myVariantId){
        variants{ id }
      }
    }
    """ 

    # Set variables object of arguments to be passed to endpoint
    variables = { "myVariantId": variant_id}

    # Set base URL of Genetics Portal GraphQL API endpoint
    base_url = "https://api.genetics.opentargets.org/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})

    # Transform API response into JSON 
    api_response_as_json = json.loads(r.text)

    # get ID
    try:
      vid = api_response_as_json['data']['search']['variants'][0]['id']
      return vid
    except:
      return None


def query_OpenTargetGenetics_genesForVariant(variant_id):
    
    # Build query string
    query_string = """
    query genesForVariant ($myVariantId: String!){
      genesForVariant(variantId:$myVariantId) {
        gene {
            symbol
            }
        overallScore
      }   
    }
    """ 

    # Set variables object of arguments to be passed to endpoint
    variables = { "myVariantId": variant_id}

    # Set base URL of Genetics Portal GraphQL API endpoint
    base_url = "https://api.genetics.opentargets.org/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})

    # Transform API response into JSON 
    api_response_as_json = json.loads(r.text)

    # get max scoring gene
    max_score = 0
    max_gene = None
    for x in api_response_as_json["data"]["genesForVariant"]:
      if x["overallScore"] > max_score:
        max_score = x["overallScore"]
        max_gene = x["gene"]["symbol"]

    return max_gene

# load files
input = sys.argv[1]
output = sys.argv[2]

# load variant list
df = pd.read_csv(input)

# map variants to gene
mapping = dict()
for index, row in df.iterrows():
  vid = query_OpenTargetGenetics_variantForRSID(row['rsid'])
  if vid is None:
    continue
  d = query_OpenTargetGenetics_genesForVariant(vid)
  if d is None:
    continue
  mapping[row['key']] = d

f = open(output, "w")
f.write("key\tgene\n")
for variant in mapping:
  f.write(variant + "\t" + mapping[variant] + "\n")
f.close()