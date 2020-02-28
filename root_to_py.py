"""Create a .npy array with the relevant data from the input .root file."""
import argparse
import numpy as np
from matplotlib import pyplot as plt


def parse_inputs():
    """Parse the command line inputs and return them in a dict, if used."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-inputFilename",
                        help="path of the input .root file with the Po data",
                        type=str)
    parser.add_argument("-minDst",
                        help="Cut minDst <= DstNumber. Default 209",
                        type=int)
    parser.add_argument("-maxDst",
                        help="Cut DstNumber <= maxDST. Default 1000",
                        type=int)
    parser.add_argument("-minCharge",
                        help="Cut minCharge <= Charge_Geo. Default 150",
                        type=int)
    parser.add_argument("-maxCharge",
                        help="Cut Charge_Geo <= maxCharge. Default 270",
                        type=int)
    parser.add_argument("-MLPCut", help="Value of the MLP cut. Default 0.3.",
                        type=float)

    outDict = vars(parser.parse_args())
    return {key: val for key, val in outDict.items() if val is not None}


def getData(tree, minDst, maxDst, minCharge, maxCharge, MLPCut):
    """Read the .root tree, set cuts on dst and MLP, return relevant vars."""
    # Define useful variables
    nEntries = tree.GetEntries()
    f_nEntries = float(nEntries)
    data = np.zeros((nEntries, 6))

    # Loop over events
    for i in xrange(nEntries):
        if i % 100000 == 0:
            print "{:.2f} percent done".format(i / f_nEntries * 100.)
        tree.GetEvent(i)
        # Apply cuts
        charge = tree.Charge_Geo
        dst = tree.DstNumber
        if dst < minDst or dst > maxDst:
            continue
        if tree.MLPv8 > MLPCut:
            continue
        if charge < minCharge or charge > maxCharge:
            continue

        # Extract relevant variables
        data[i, :] = tree.x, tree.y, tree.z, charge, tree.RealTime, dst

    # Reshape data to remove empty elements
    # TODO: Could maybe replace by filling lists, then converting to np array?
    m = data[:, 0] != 0
    newData = np.zeros((np.sum(np.array(m, dtype=int)), 6))
    for i in xrange(data.shape[1]):
        newData[:, i] = data[:, i][m]

    return newData


def main(inputFilename="Po_201601-20200902_aligned_1902.root", minDst=209,
         maxDst=1000, minCharge=150, maxCharge=270, MLPCut=0.3):
    """Read the tree, then save the extracted data to .npy array."""
    # Read in the tree
    f = TFile(inputFilename)
    tree = f.Get("t")

    # Get the data
    data = getData(tree, minDst, maxDst, minCharge, maxCharge, MLPCut)
    np.save("data_dst_{}_{}_mlp_{}.npy".format(minDst, maxDst, minCharge, maxCharge, MLPCut), data)
    print data.shape


if __name__ == '__main__':
    args = parse_inputs()
    # Have to import ROOT here, otherwise it messes with argparse...
    from ROOT import TFile
    main(**args)
