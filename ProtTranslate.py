from Bio.Seq import translate
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio import Seq
import shutil

import logomaker
import pandas as pd
import matplotlib.pyplot as plt

from ambigMap import ambigMap

# def generateExternalMap():
    # # import json

    # # with open("ambigMap.json", "w") as f:
    # #     json.dump(ambigMap, f, indent=4, separators=(",", ": "))

    # with open("ambigMap1.py", "w") as f:
    #     f.write("ambigMap = ")
    #     f.write(repr(ambigMap))
        
    # import pprint

    # with open("ambigMap.py", "w") as f:
    #     f.write("ambigMap = \\\n")
    #     pprint.pprint(ambigMap, stream=f, indent=2, width=80)

def textBool(bool):
    return "Yes" if bool else "No"

def parseInput(fileName):
    with open(fileName, 'r') as inputFile:
        lines = [line.strip() for line in inputFile if line.strip()]
    return str(lines[0]), str(lines[1])
    
def findIndex(ints, best):
    for i in range(len(ints)):
        if ints[i] == best:
            return i+1 
    return 0
    
def translatorStrat(stratNum, codingDNA,protStrand, endStops, stopSymbol, looseHam=False):
    options = []
    scores = []
    keepStrat=False
    max_table_id = max(CodonTable.generic_by_id.keys())
    for i in range(1, max_table_id + 1):
        translStr = tryTranslate(seq=codingDNA, t=i, stopSymbol=stopSymbol, endStops=endStops)
        if translStr is None:
            continue
        
        print(f"Table {i}: ",end="")
        if len(translStr) <= 81:
            print(f"{translStr}")
            
        hamm= looseHamDist(translStr,protStrand) if looseHam else hammingDist(translStr,protStrand)
        print("hammingDist: ", hamm)
       
        if protStrand == translStr:
            options.append(i)
            keepStrat=True
            
        scores.append(hamm)  
        printDivider('=')
    if not scores:
        print(f"[!] Strategy {stratNum} failed: no valid translations.")
        return
    minScoreIndex = findIndex(scores, min(scores))
    
    print(f"Strat {stratNum}:")
    print(min(options)) if options else (print(minScoreIndex) if scores else print("error: No options"))
    printDivider('=')
    printDivider('=')
    return {
        "strategy": stratNum,
        "best_table": min(options) if options else minScoreIndex,
        "min_score": min(scores),
        "matched_exactly": keepStrat,
    }
    
    # if keepStrat:
    #     return
    
def strDivider(char='='):
    try:
        width = shutil.get_terminal_size().columns
    except:
        width = 80  # fallback if anything fails
    return(char * width)    
    
def printDivider(char='='):
    print(strDivider(char))
    
def printToFile(string):
    fileName = "probeOut.txt"
    with open(fileName, 'a') as of: # a for append instead of w for write
        of.write(string+'\n')
    

def findCodonTable(codingDNA, protStrand):
    #Strat 1
    translatorStrat(1, codingDNA, protStrand, endStops=False, stopSymbol='', looseHam=False)
    
    #Strat 2
    translatorStrat(2, codingDNA, protStrand, endStops=False, stopSymbol='*', looseHam=False)
    
    #Strat 3
    translatorStrat(3, codingDNA, protStrand, endStops=False, stopSymbol='*', looseHam=True)
    
    #Strat 4
    translatorStrat(4, codingDNA, protStrand, endStops=True, stopSymbol='*', looseHam=True)
    
    
def probeTables(codons = ["TAG","TAA","TGA"],stopSymbol=''):
    open("probeOut.txt", 'w').close()  # Clear file at the start
    for codon in codons:
        print(f"For Codon: {codon}")
        printToFile(f"For Codon: {codon}")
        # for i in range(1,7):
        max_table_id = max(CodonTable.generic_by_id.keys())
        for i in range(1, max_table_id + 1):
            # printToFile(f"{i}")
            if i == 7 or i == 8 or (17 <= i and i <= 20): #skip Broken Tables in tool
                continue
            # translStr = translate(codon,table=i, to_stop=False, stop_symbol=stopSymbol)
            translStr = tryTranslate(seq=codon, t=i, stopSymbol=stopSymbol, endStops=False, toFile=True)
            if translStr is None:
                continue
            print(f"Table {i}: {translStr}")
            printToFile(f"Table {i}: {translStr}")
        printDivider()
        printToFile(strDivider())

def tryTranslate(seq, t, stopSymbol, endStops=False, toFile=False):
    try:
        return str(translate(seq, table=t, stop_symbol=stopSymbol, to_stop=endStops))
    except Exception as e:
        print(f"[!] Table {t} caused error: {e}") if not toFile else printToFile(f"[!] Table {t} caused error: {e}")
        return None

def hammingDist(A,B):
    differences = []
    d= 0 if len(A) == len(B) else abs(len(A) - len(B))
    if len(A) != len(B):
        print(f"Different Lengths {len(A)} != {len(B)} dif = {abs(len(A) - len(B))}")   
    for i in range(min(len(A),len(B))):
        if A[i] != B[i]:
            d+=1
            # print
            # differences.append(f"{A[i]} != {B[i]}")
    return d
    # return d if len(differences) == 0 else (d, differences)

def looseHamDist(A,B):
    d=0
    if len(A) != len(B):
        print(f"Different Lengths {len(A)} != {len(B)}, dif = {abs(len(A) - len(B))}")
    for i in range(min(len(A),len(B))):
        if A[i] != B[i] and not(
            (A[i] == '*' and B[i] in ['Q', 'W']) or 
            (B[i] == '*' and A[i] in ['Q', 'W'])
            ):
            d+=1  
    return d 

def stitch(strs):
    newStr= ""
    for string in strs:
        newStr+=string
    return newStr

def generateCodonStr():
    return stitch(generateCombos(3))

def generateCombos(n):
    if n == 1:
        return ['A','C','G','T']
    else:
        newCombos = []
        for combo in generateCombos(n-1):
            for letter in ['A','C','G','T']:
                newCombos.append(combo+letter)
        return newCombos
            
# if "[!]" in  outputStr or  "error" in  outputStr: #skip tables 7-8, 17-20
#     continue
def parseOutput(outputBlock):
    codonMap = dict()
    aminoMap = dict()
    codon = None
    for outputStr in outputBlock:
        outputStr = outputStr.strip()
        if "Codon" in outputStr:
            codon = outputStr.split(":")[-1].strip()
            aminoMap = {}
            if codon not in codonMap:
                codonMap[codon] = aminoMap
        elif "Table" in outputStr and codon:
            amino = outputStr[-1]
            tableNum = int(outputStr.split(" ")[1][:-1])
            if amino not in aminoMap:
                aminoMap[amino] = [tableNum]
            else:
                aminoMap[amino].append(tableNum)
    if codon:
        codonMap[codon] = aminoMap
    return codonMap
        
def outputParser():
    outputBlocks = []
    outputBlock = []
    with open("probeOut.txt", 'r') as inputFile: 
        lines = inputFile.readlines()
        for line in lines:
            if '=' in line:
                outputBlocks.append(outputBlock)            
                outputBlock = []
            else:
                outputBlock.append(line)
    return outputBlocks
        
def checkAllCodons():
    codonMap = dict()
    for outputBlock in outputParser():
        blockResult = parseOutput(outputBlock)
        for codon, aminoMap in blockResult.items():
            if codon not in codonMap:
                codonMap[codon] = dict()
            for amino, tables in aminoMap.items():
                if amino not in codonMap[codon]:
                    codonMap[codon][amino] = tables.copy()
                else:
                    codonMap[codon][amino].extend(tables)
                    # codonMap[codon][amino][1] += count
    return codonMap

def checkForAmb():
    codonMap = checkAllCodons()
    ambigCodons=dict()
    for codon, aminoMap in codonMap.items():
        if len(aminoMap) > 1:
            ambigCodons[codon] = aminoMap  # keep the entire aminoMap
    return ambigCodons

def printAmbigMap(ambigCodons):
    for codon, aminoMap in ambigCodons.items():   
        print(f"{codon}: ")
        for amino, locations in aminoMap.items():
            freq = len(locations)
            print(f"{amino}: {freq} times in ", end="")
            for tableNum in locations:
                print(f"{tableNum} ", end="")
            print()
                   
def genAmbigList(ambigCodons):
    codonList = []
    for codon, aminoMap in ambigCodons.items():
        codonList.append(codon)
        
    return codonList

def extractFreqs(ambigMap, Codon):
    return{amino: [len(locations)] for amino, (locations) in ambigMap[Codon].items()}

def logoMaker(Codon):
    freqMap = extractFreqs(ambigMap, Codon)
    allAminos = "ACDEFGHIKLMNPQRSTVWY"
    
    dataFrame = pd.DataFrame(freqMap)

    # Fill in all 20 amino acids so logomaker doesn’t break
    for amino in 'ACDEFGHIKLMNPQRSTVWY*':
        if amino not in dataFrame.columns:
            dataFrame[amino] = 0

    # Reorder columns
    dataFrame = dataFrame[['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*']]

    # Normalize to get frequency
    dataFrame = dataFrame.div(dataFrame.sum(axis=1), axis=0)

    # Make the logo
    logomaker.Logo(dataFrame, color_scheme='chemistry', font_name='Arial', shade_below=.5, fade_below=.5)
    plt.title(f"Codon {Codon} - Amino Acid Mapping Profile")
    plt.show()

ambigCodons = ["AAA","AGA","AGG","ATA","CTA","CTC","CTG","CTT","TAA","TAG","TCA","TGA","TTA"]
def isAmbig(codon):
    return codon in ambigCodons

def hasFrameShift(DNA,Aminos):
    print(f"DNA:    {(len(DNA) / 3)}")
    print(f"Aminos: {(len(Aminos))}")
    return (len(DNA) / 3) != (len(Aminos)) 

def printBool(bool):
    print(textBool(bool))

#templates
    # aminoMap = { 'A': [[tableNums], freq]} #old form
    # aminoMap = { 'A': [tableNums]}        #refactor form
    # codonMap = {"AAA" : [aminoMap]}

# main# main# main# main# main# main# main# main# main# main# main# main# main# main# main
# inputFile = "rosalind_ptra.txt"
# #"translateSample.txt" 

# Sample Set
DNA = "ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
Aminos = "MAMAPRTEINSTRING"
# Sample Output
# 1

# DNA,Aminos = parseInput(inputFile)

# findCodonTable(DNA,Aminos)

# print(generateCodonStr())

# # probeTables()
# probeTables(generateCombos(3),'*')


# print(checkAllCodons())
ambigMap = checkForAmb()
printAmbigMap(ambigMap)
print(genAmbigList(ambigMap))

printBool(hasFrameShift(DNA,Aminos))

logoMaker(Codon="AAA")

