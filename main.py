# script to extract information about extracellular
# tyrosine residues from the Uniprot database

import json, sys, requests, csv

# function to return a dictionary response
def getData(accession):
    requestURL = "https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession=" + accession + "&types=TOPO_DOM%2CCHAIN"
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
    i = 0 # count number of extracellular tyrosine residues
    j = 0 # count number of extracellular amino acids
    if len(topoArr) == 1:
        topoDict = topoArr[0]
        roi = sequence[int(topoDict['begin'])-1:int(topoDict['end'])]
        i += roi.count('Y')
        j += len(roi)
    else:
        for topoDict in topoArr:
            if (topoDict['description'] in features):
                roi = sequence[int(topoDict['begin'])-1:int(topoDict['end'])]
                i += roi.count('Y') # count number of extracellular tyrosine residues
                j += len(roi) # count number of extracellular amino acids
    return [i,j]

# get file arguments
files = sys.argv
accessionFile = files[0]
dataFile = files[1]

# get list of accession values
with open(accessionFile, newline='') as f:
    reader = csv.reader(f)
    accessionList = list(reader)[0]

accessionList = ["O88944"]

features = "Extracellular"
# dict_keys(['accession', 'entryName', 'sequence', 'sequenceChecksum', 'taxid', 'features'])

tyrData = []

for accession in accessionList:
    responseDict = getData(accession)
    sequence = responseDict['sequence'] # the peptide sequence
    totalTyr = sequence.count('Y') # the total number of tyrosines
    topoArr = responseDict['features'] # an array with the topological features as dicts
    extracellularTyr = countTyr(sequence, topoArr, features)
    description = topoArr[0]['description']
    tyrData.append([accession,description,extracellularTyr[0],extracellularTyr[1],len(sequence)])

# write data to file
with open(dataFile, 'wt') as fp:
    writer = csv.writer(fp, delimiter=',')
    for i in tyrData:
        writer.writerow(i)
