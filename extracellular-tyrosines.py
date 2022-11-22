# script to extract information about extracellular
# tyrosine residues from the Uniprot database

import json, sys, requests, csv, numpy as np

# function to return a dictionary response from the Uniprot database given an accession number
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

# function to read in a numpy matrix from csv
def csvToMat(filename):
    with open(filename, newline='') as f:
        reader = csv.reader(f)
        mat = np.array(list(reader))
    return mat

# function to write a matrix to a csv
def csvFromMat(outputFile, outputMat):
    with open(outputFile, 'wt') as fp:
        writer = csv.writer(fp, delimiter=',')
        for i in outputMat:
            writer.writerow(i)

# get file arguments
files = sys.argv
topHitsFile = files[1] # top hits file path
PSMFile = files[2] # PSM data file path
#dataFile = files[3] # output file path
features = "Extracellular"

accessionMat = csvToMat(topHitsFile) # read in top hits
PSMMat = csvToMat(PSMFile) # read in PSM data

uniqueAcc = list(set([j for i in range(0,5) for j in accessionMat[2:,(0+3*i)] if j]))

# uniqueAcc = uniqueAcc.tolist()
# #uniqueAcc.remove("Q3UH53") # unknown positions in database
# #uniqueAcc.remove("Q8BLK3") # unknown positions in database

# get all tyrosine and extracellular residue data
uniqueAccData = []
# ["Accession", "number of extracellular tyrosines", "number of extracellular amino acids","total length of amino acid"]
for accession in uniqueAcc:
    responseDict = getData(accession)
    sequence = responseDict['sequence'] # the peptide sequence
    totalTyr = sequence.count('Y') # the total number of tyrosines
    topoArr = responseDict['features'] # an array with the topological features as dicts
    extracellularTyr = countTyr(sequence, topoArr, features)
    description = topoArr[0]['description']
    uniqueAccData.append([accession,extracellularTyr[0],extracellularTyr[1],len(sequence)])

# put unique accession numbers and their tyrosine data in a dictionary
uniqueAccDict = dict()
for i in uniqueAccData:
    uniqueAccDict[i[0]] = i[1:]

# put all PSM data in a dictionary keyed by accession numbers
PSMDict = dict()
for i in PSMMat[4:]:
    if i[0] in uniqueAcc:
        div4arr = [j for j in [i[5],i[7],i[9]]]
        div7arr = [j for j in [i[11],i[13],i[15]]]
        div14arr = [j for j in [i[17],i[19],i[21]]]
        div21arr = [j for j in [i[23],i[25],i[27]]]
        div28arr = [j for j in [i[29],i[31],i[33]]]
        PSMDict[i[0]] = dict(div4=div4arr,div7=div7arr,div14=div14arr,div21=div21arr,div28=div28arr)

divs = ["div4","div7","div14","div21","div28"]
divDataDict = dict() # the dictionary with all data for output
for i in divs:
    divDataDict[i] = dict(accArr=[],tyrMat=[],PSMArr=[],names=[])

for i in range(0,5):
    divDataDict[divs[i]]['accArr'] = [j for j in accessionMat[2:,(0 + 3*i)] if j]
    divDataDict[divs[i]]['names'] = [j for j in accessionMat[2:,(1 + 3*i):(3 + 3*i)] if j[0]]
    divDataDict[divs[i]]['tyrMat'] = [uniqueAccDict[j] for j in divDataDict[divs[i]]['accArr']]
    divDataDict[divs[i]]['PSMArr'] = [np.mean([float(k) for k in PSMDict[j][divs[i]]]) for j in divDataDict[divs[i]]['accArr']]

# for each accession: accession number, protein name, gene, num extracellular tyrosines, num extracellular amino acids, total length of amino acids, PSM
header = ["acc","name","gene","#tyr","#ExC AA","#AA","PSM"]

outputDict = dict()
for i in divs:
    outputDict[i] = np.concatenate(([divDataDict[i]['accArr']], np.array(divDataDict[i]['names']).transpose(), np.array(divDataDict[i]['tyrMat']).transpose(), [divDataDict[i]['PSMArr']])).transpose()
    # write output to file
    outputFileName = i + "outputFile.csv"
    csvFromMat(outputFileName,np.concatenate(([header],outputDict[i])))
