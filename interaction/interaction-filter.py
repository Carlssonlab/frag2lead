#!usr/bin/env python
__author__  = 'Andreas Johan Luttens'
__date__    = '2019-12-31'
__email__   = 'andreas.luttens@gmail.com'
__version__ = '0.1'

"""
Extract poses for compounds that engage in a specific interaction with parsed residues.

Caveat emptor! Script assumes single pose per compound!

[Andreas Luttens, 2024]
"""

import argparse
import os
import gzip
from openeye import oechem

def ParseArgs():
    parser = argparse.ArgumentParser(epilog='https://github.com/carlssonlab/frag2lead')
    parser.add_argument('-i', '--input_file', 
                        help = 'coordinate file input molecules', 
                        required=True)
    parser.add_argument('-o', '--output_file', 
                        help = 'coordinate file output molecules', 
                        required=True)
    parser.add_argument('-p', '--protein_file', 
                        help = 'coordinate file proteine', 
                        required=True)
    parser.add_argument('-s', '--interactionString', 
                        help = 'logic string containing required interaction patterns', 
                        required=True)
    args = parser.parse_args()
    return args

def ProcessSubString(string):
    residueNumber = ""
    typeString = ""

    for c in string:
        if c.isnumeric():
            residueNumber += c
        elif c.isalpha:
            typeString += c

    return int(residueNumber), typeString

def ProcessInteractionString(interactionString):
    """
    Process string containing specific interactions, return residues and interaction predicate.
    """ 
    
    # OEHBondInteractionHintType_LigandDonates

    # library with possible interactions
    interactionLibrary = {"B":oechem.OEIsSaltBridgeInteractionHint(), 
                          "H":oechem.OEIsIntermolecularHBondInteractionHint(), 
                          "HD":oechem.OEHasInteractionHintType(oechem.OEHBondInteractionHint(oechem.OEHBondInteractionHintType_LigandDonates)),
                          "HA":oechem.OEHasInteractionHintType(oechem.OEHBondInteractionHint(oechem.OEHBondInteractionHintType_LigandAccepts)),
                          "S":oechem.OEIsStackingInteractionHint(), 
                          "C":oechem.OEIsCationPiInteractionHint(),
                          "X":oechem.OEIsHalogenBondInteractionHint()
                          }

    # get the specific interaction predicates
    predicateLibrary = {}
    orVector = []
    orAll = []

    rawInteractions = interactionString.split("-")
    for interaction in rawInteractions:
        if ":" in interaction:
            # add an empty list to store residue numbers
            orVector.append([])
            for option in interaction.split(":"):
                residueNumber, typeString = ProcessSubString(option)
                # if the residue is already required to have a specific interaction, merge the predicates
                if residueNumber in predicateLibrary.keys():
                    predicateLibrary[residueNumber] = oechem.OEAndInteractionHint(predicateLibrary[residueNumber], interactionLibrary[typeString])
                else:
                    predicateLibrary[residueNumber] = interactionLibrary[typeString]
                orVector[-1].append(residueNumber)
                orAll.append(residueNumber)
        else:
            residueNumber, typeString = ProcessSubString(interaction)

            # if the residue is already required to have a specific interaction, merge the predicates
            if residueNumber in predicateLibrary.keys():
                predicateLibrary[residueNumber] = oechem.OEAndInteractionHint(predicateLibrary[residueNumber], interactionLibrary[typeString])
            else:
                predicateLibrary[residueNumber] = interactionLibrary[typeString]

    residues = predicateLibrary.keys()

    return residues, predicateLibrary, orVector, orAll

def IsInteractor(proteinMol, candidateMol, residues, predicateLibrary, orVector, orAll):
    """
    Perceive intermolecular interactions, check if predicates are met.
    """

    # create bindingsite object to perceive interactions
    asite = oechem.OEInteractionHintContainer(proteinMol, candidateMol)
    if not oechem.OEIsValidActiveSite(asite):
        oechem.OEThrow.Fatal("Cannot initialize interaction site!")
    oechem.OEPerceiveInteractionHints(asite)
    
    # initiate interaction library
    interactionLibrary= {}

    # iterate over residues, that engage in an interaction
    for res in oechem.OEGetResidues(asite.GetMolecule(oechem.OEProteinInteractionHintComponent())):
        # if the residue is in the list of interest
        residueNumber = res.GetResidueNumber()
        if residueNumber in residues:
            # check if interaction is specific AND with this residue
            if sum([1 for i in asite.GetInteractions(oechem.OEAndInteractionHint(oechem.OEHasResidueInteractionHint(res), predicateLibrary[residueNumber]))]) > 0:
                interactionLibrary[residueNumber] = 1
            else:
                interactionLibrary[residueNumber] = 0
        # skip residue
        else:
            pass
    
    # check if all requirements are met
    # check if at least one interaction is found in all or-statements
    if 0 in [sum([interactionLibrary[residueNumber] for residueNumber in orVector[index]]) for index,couple in enumerate(orVector)]:
        return False
    else:
        if 0 in [interactionLibrary[residueNumber] for residueNumber in residues if residueNumber not in orAll]:
            return False
        else:
            return True

def CurateViewDockHeaders(inputFile, outputFile):
    """
    Change order of header and block of mol2 text
    """

    bufferString = ""
    status = False
    with open(inputFile, 'r') as inhandle:
        with open(outputFile, 'w') as outhandle:
            for line in inhandle.readlines():
                if line.startswith("@<TRIPOS>MOLECULE"):
                    status = True
                elif line.startswith("##########"):
                    status = False
                elif line.startswith("@<TRIPOS>ATOM"):
                    outhandle.write(bufferString)
                    bufferString = ""
                if status:
                    bufferString += line
                else:
                    outhandle.write(line)
        outhandle.close()
    inhandle.close()

def RetrieveHeaders(inputFile, candidates):
    """
    For a given inputfile, grab the headers in the mol2 file and store them in an index dictionary.
    """

    headerDictionary = {}
    energyDictionary = {}
    countDictionary = {}
    instanceDictionary = {}
   
    header = ""
    status = False
    
    # if all candidates are found, exit, tricky part
    count=0
    total=len(candidates)
    seen = set()

    # initialize in the candidate counters
    for name in candidates:
        countDictionary[name] = 0
    
    currentName=None
    readStatus = True
    # read binary test.mol2.gz
    with gzip.open(inputFile, 'rt') as inhandle:
        # don't keep all lines in memory
        # LSD filed will be very big
        while readStatus:
            line = inhandle.readline()
            if line:
                if line.startswith("##########"):
                    status = True
                    header += line
                    if "  Name" in line:
                        name = line.split()[-1].strip()
                        # if all candidates have been seen and new name
                        if len(seen)==total and name!=currentName:
                            inhandle.close()
                            return headerDictionary, instanceDictionary
                        currentName=name
                    elif "Total Energy" in line:
                        value = float(line.strip().split()[-1])
                else:
                    if status:
                        if name in candidates:
                            countDictionary[name]+=1
                            if name in headerDictionary:
                                if energyDictionary[name] > value:                            
                                    headerDictionary[name] = header
                                    energyDictionary[name] = value
                                    instanceDictionary[name]=countDictionary[name]
                            else:
                                headerDictionary[name] = header
                                energyDictionary[name] = value
                                instanceDictionary[name]=countDictionary[name]

                            seen.add(name)

                        header = "" 
                        status = False
            else:
                readStatus = False
            
        inhandle.close()
    return headerDictionary, instanceDictionary

def main():
    """
    Driver function
    """

    # parse the arguments
    args = ParseArgs()
    
    # read the protein structure
    proteinMol = oechem.OEGraphMol()
    proteinStream = oechem.oemolistream()
    if not proteinStream.open(args.protein_file):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % args.protein_file)
    if not oechem.OEReadMolecule(proteinStream, proteinMol):
        oechem.OEThrow.Fatal("Unable to read molecule from %s" % args.protein_file)
    proteinStream.close()
   
    # make sure residues are read in
    if not oechem.OEHasResidues(proteinMol):
        oechem.OEPerceiveResidues(proteinMol, oechem.OEPreserveResInfo_All)

    # process parsed interaction string
    residues, predicateLibrary, orVector, orAll = ProcessInteractionString(args.interactionString)

    # read the inputfile as a stream
    inputfileStream = oechem.oemolistream(args.input_file)
    if not inputfileStream.open(args.input_file):
        oechem.OEThrow.Fatal("Unable to open %s inputfile for reading" % args.input_file)

    # file to store output names
    outputStream = open(args.output_file, 'w')
    
    # iterate over all molcules in the inputfile   
    for molecule in inputfileStream.GetOEGraphMols():
        molname = molecule.GetTitle().split()[0]
        
        if IsInteractor(proteinMol, molecule, residues, predicateLibrary, orVector, orAll):
            outputStream.write(f'{molname}\n')

    outputStream.close()
    
if __name__ == '__main__':
    main()

