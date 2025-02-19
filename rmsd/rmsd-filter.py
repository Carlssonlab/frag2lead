#!usr/bin/env python

__author__  = 'Andreas Johan Luttens'
__date__    = '2019-12-31'
__email__   = 'andreas.luttens@gmail.com'
__version__ = '0.1'

"""
Extract poses from docking screen within RMSD cutoff with crystal structure.
RMSD is based on heavy atom common core calculations.

[Andreas Luttens, 2024]
"""

import os
import argparse
from openeye import oechem
from math import sqrt
import gzip

def ParseArgs():
    """
    Parse the arguments to the main driver function
    """
    parser = argparse.ArgumentParser(epilog='https://github.com/carlssonlab/frag2lead')
    parser.add_argument('-r', '--reference', 
                        help = 'file containing reference molecule', 
                        required=True)
    parser.add_argument('-c', '--comparison', 
                        help = 'file containg poses to reference', 
                        required=True)
    parser.add_argument('-s', '--smarts', 
                        help = "SMARTS pattern to match", 
                        required=True)
    parser.add_argument('-t', '--threshold', 
                        help = 'RMSD threshold',
                        type = float, 
                        default = 2.0)
    parser.add_argument('-o', '--outputfile', 
                        help = 'file where output poses will be stored',
                        required=True)
    args=parser.parse_args()
    return(args)

def breakdownHeaders(inputfile):
    """
    For a given inputfile, grab the headers in the mol2 file and store them in an index dictionary.
    """

    headerDictionary = {}
    energyDictionary = {}

    index = 0
    header = ""
    status = False

    with gzip.open(inputfile, 'rb') as filehandle:
        for line in filehandle.readlines():
            line = line.decode("utf-8")
            if line.startswith("##########"):
                status = True
                header += line
                if "  Name" in line:
                    name = line.split()[-1].strip()
                    if name not in energyDictionary:
                        energyDictionary.setdefault(name, {})
                elif "Total Energy" in line:
                    energy = float(line.split()[-1].strip())
                    energyDictionary[name][index]=energy
            else:
                if status:
                    headerDictionary[index] = header
                    header = ""
                    index += 1
                    status = False
        filehandle.close()

    return headerDictionary, energyDictionary

def calculateRMSD(refmol, fitmol, overlay):
    """
    Simple RMSD calculation based on common atom predicates.
    """
    atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Aromaticity
    bondexpr = 0
    mcss = oechem.OEMCSSearch(oechem.OEMCSType_Exhaustive)
    mcss.Init(refmol, atomexpr, bondexpr)
    mcss.SetMCSFunc(oechem.OEMCSMaxBondsCompleteCycles())

    unique = True
    overlay = False
    for match in mcss.Match(fitmol, unique):
        rmsd =  oechem.OERMSD(mcss.GetPattern(), fitmol, match, overlay)
        return rmsd

def reorderFile(inputfile, outputfile):
    """
    Change order of header and block of mol2 text
    """

    buffer = ""
    status = False
    with open(inputfile, 'r') as inhandle:
        with open(outputfile, 'w') as outhandle:
            for line in inhandle.readlines():
                if line.startswith("@<TRIPOS>MOLECULE"):
                    status = True
                elif line.startswith("##########"):
                    status = False
                elif line.startswith("@<TRIPOS>ATOM"):
                    outhandle.write(buffer)
                    buffer = ""
                if status:
                    buffer += line
                else:
                    outhandle.write(line)
        outhandle.close()
    inhandle.close()

def main():
    """
    Driver function
    """

    # parse the arguments
    args = ParseArgs()

    # read the reference molecule
    refmol = oechem.OEMol()
    referenceStream = oechem.oemolistream(args.reference)
    if not referenceStream.open(args.reference):
        oechem.OEThrow.Fatal("Cannot open %s input file!" % args.reference)
    oechem.OEReadMolecule(referenceStream, refmol)
    oechem.OEAddExplicitHydrogens(refmol)
    referenceStream.close()

    # testing
    comparisonStream = oechem.oemolistream(args.comparison)
    if not comparisonStream.open(args.comparison):
        oechem.OEThrow.Fatal("Cannot open %s input file!" % args.comparison)

    # create outputfile
    outputStream = oechem.oemolostream("intermediate.mol2")
    if not outputStream.open("intermediate.mol2"):
        oechem.OEThrow.Fatal("Cannot open %s output file!" % "intermediate.mol2")

    # create a substructure search object
    ss = oechem.OESubSearch(args.smarts)  
    ss.SetMaxMatches(1)

    oechem.OEPrepareSearch(refmol, ss)
    refatoms = None

    refcoords = refmol.GetCoords()
    subrefcoords = None
    for idx,match in enumerate(ss.Match(refmol,True)):
        if idx>0:
            break
        else:
            subrefcoords = [refcoords[a.GetIdx()] for a in match.GetTargetAtoms()]

    # obtain headers for every molecular object
    headerDictionary, energyDictionary = breakdownHeaders(args.comparison)

    # store indices of valid compound poses
    validPoses = {}

    # store indices of best compound poses
    bestPoses = []

    # iterate over molecular objects
    objectCounter = 0
    currentName = ""
    bestEnergy = 9999.99
    bestPose = ""

    for molobj in comparisonStream.GetOEGraphMols():
        # get the match between reference and this molecule
        oechem.OEPrepareSearch(molobj, ss)

        rmsdobj = 0.0
        natoms = 0
        submolcoords = None
        for idx,match in enumerate(ss.Match(molobj,True)):
            if idx > 0:
                break
            else:
                molcoords = molobj.GetCoords()
                submolcoords = [molcoords[a.GetIdx()] for a in match.GetTargetAtoms()]
        
        if submolcoords is None:
            print("SMARTS did not match: %s %s"%(oechem.OEMolToSmiles(molobj), molobj.GetTitle()))
            continue

        for atom in range(len(subrefcoords)):
            # x y z
            for i in range(3):
                rmsdobj += (subrefcoords[atom][i] - submolcoords[atom][i])**2
            natoms += 1
        rmsdobj /= natoms
        rmsdobj = sqrt(rmsdobj)
        
        if rmsdobj > 0.0 and rmsdobj <= args.threshold:
            objname = molobj.GetTitle().split()[0]
            objenergy = energyDictionary[objname][objectCounter]
            if objname != currentName:
                bestPoses.append(bestPose)
                currentName = objname
                bestEnergy = objenergy
                bestPose = objectCounter

            else:
                if objenergy < bestEnergy:
                    bestEnergy = objenergy
                    bestPose = objectCounter

        # always update counter
        objectCounter +=1
    
    comparisonStream.close()

    comparisonStream = oechem.oemolistream(args.comparison)
    if not comparisonStream.open(args.comparison):
        oechem.OEThrow.Fatal("Cannot open %s input file!" % args.comparison)

    bestPoses.append(bestPose)

    objectCounter = 0
    for molobj in comparisonStream.GetOEGraphMols():
        if objectCounter in bestPoses:
            oechem.OESetComment(molobj, headerDictionary[objectCounter])
            oechem.OEWriteMolecule(outputStream, molobj)
        objectCounter += 1

    comparisonStream.close()
    outputStream.close()

    reorderFile("intermediate.mol2", args.outputfile)

    os.remove("intermediate.mol2")

if __name__ == '__main__':
    main()
