# script to extract information about extracellular
# tyrosine residues from the Uniprot database

import json, sys, requests

# function to return a dictionary response
def getData(accession):
    requestURL = "https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession=" + accession + "&types=TOPO_DOM"
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = r.text
    responseDict = json.loads(responseBody)[0]
    return responseDict

# function to count tyrosines in a sequence given
# an array of topological features and an array
# of features of interest
def countTyr(sequence, topoArr, features):
    i = 0
    for topoDict in topoArr:
        if (topoDict['description'] in features):
            # print(topoDict['description'] + " is in " + features)
            # print("Finding extracellular tyrosines from position " + topoDict['begin'] + " to position " + topoDict['end'])
            # print("The sequence of interest is:")
            # print(sequence[int(topoDict['begin'])-1:int(topoDict['end'])-1])
            i += sequence[int(topoDict['begin'])-1:int(topoDict['end'])-1].count('Y')
    return i

accession = "Q810U4" #"Q3UHK6" #"Q810U3"
features = "Extracellular"
# dict_keys(['accession', 'entryName', 'sequence', 'sequenceChecksum', 'taxid', 'features'])
responseDict = getData(accession)
sequence = responseDict['sequence'] # the peptide sequence
totalTyr = sequence.count('Y') # the total number of tyrosines
topoArr = responseDict['features'] # an array with the topological features as dicts

extracellularTyr = countTyr(sequence, topoArr, features)
print("Extracellular tyrosine residues: " f'{extracellularTyr}')
print("Total tyrosine residues: " + f'{totalTyr}')
print(responseDict)
