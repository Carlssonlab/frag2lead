#!/usr/bin/env python

__author__  = 'Andreas Johan Luttens'
__date__    = '2019-12-05'
__email__   = 'andreas.luttens@gmail.com'
__version__ = '0.1'

"""
Retrieve compounds that match your SMARTS query in the Enamine REAL set.
 
[Andreas Luttens, 2024]
"""

import argparse
from openeye.oechem import *
from threading import Thread

def parseArgs():
    """
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(epilog='https://github.com/carlssonlab/frag2lead')
    parser.add_argument('-i', '--input', 
                        help = "set of molecules you wish to search for substructures", 
                        required=True)
    parser.add_argument('-p', '--pattern', 
                        help = "SMARTS to match", 
                        required=True)
    parser.add_argument('-w', '--weight', 
                        help ="maximum weight of compound", 
                        required=False, 
                        default="None")
    parser.add_argument('-r', '--reverse', 
                        help = "apply reverse logic or not", 
                        action = 'store_true', default=False)
    parser.add_argument('-o', '--output', help = "name of the output", 
                        required=True)
    parser.add_argument('-n', '--nthreads', 
                        help = "number of threads",
                        default=1, 
                        type=int) 
    parser.parse_args()
    parser.set_defaults(verbose=False)
    args = parser.parse_args()
    return args

def prepQuery(pattern, input, output):
    """
    Prepare the search engine, open streams.
    """

    # init stream targets
    inputStream = oemolithread(input)
    inputStream.SetFormat(OEFormat_CSV)
    if not inputStream.open(input):
        OEThrow.Fatal("Unable to open %s targets for reading" % input)

    # init stream output
    outputStream = oemolothread(output)
    outputStream.SetFormat(OEFormat_CSV)
    if not outputStream.open(output):
        OEThrow.Fatal("Unable to open %s output for writing" % output)

    searchEngine = OESubSearch()
    if not searchEngine.Init(pattern):
        OEThrow.Fatal("Unable to parse smarts: %s" % pattern)

    # return work objects
    return inputStream, outputStream, searchEngine

def ThreadedSearch(pattern, weight, input, output, reverse, nthreads):
    """
    Divide the search query over N threads.
    """

    # prepare the query by reading all options and settings from the MDL query file
    inputStream, outputStream, searchEngine = prepQuery(pattern, input, output)

    # Threading section
    threads = []
    threadedDots = OEThreadedDots(50000, 5000, "matches")

    class MultiCoreSearch(Thread):
        def run(self):
            # iterate over all target molecules
            for candidate in inputStream.GetOEGraphMols():
    		
                # add explicit hydrogens, might be necessary for your query
                OEAddExplicitHydrogens(candidate)
		
		            # if a weight limit is installed, do it before expensive SMARTS matching
		            if weight is not "None":
		                if OECalculateMolecularWeight(candidate) > weight:
			                  continue

                # start comparing
                OEPrepareSearch(candidate, searchEngine)

                # if there is a match or not, write out depending on the reverse flag
                if searchEngine.SingleMatch(candidate) != reverse:
                    OEWriteMolecule(outputStream, candidate)
                    threadedDots.Update()

    # Get the number of processor on current machine
    for index in range(nthreads):
        print("Started thread number %d"%(index+1))
        # split the workload over N threads
        thread = MultiCoreSearch()
        thread.start()
        threads.append(thread)

    # gather results
    for thread in threads:
        thread.join()

    threadedDots.Total()

    inputStream.close()
    outputStream.close()


def main():
    """
    Run the main script.
    """

    # retrieve the argument vector
    args = parseArgs()

    # perform a threaded substructure search
    ThreadedSearch(args.pattern, args.weight, args.input, args.output, args.reverse, args.nthreads)

if __name__ == "__main__":
    main()
