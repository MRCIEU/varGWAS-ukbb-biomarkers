import json

with open('data/opentargets/ot.json') as f:
    data = json.load(f)
    print(data)