#!/Users/martinmccullagh/anaconda/bin/python


import scipy
import sys
import getopt
import os
import numpy
import math
import MDAnalysis
import numpy.linalg
from MDAnalysis.analysis.align import *

# define some parameters
thresh = 1E-5
maxIter = 100

# read in command line arguments
inPsfFile = sys.argv[1]
coordDcdFile = sys.argv[2]

#print out log info
print "PSF file:", inPsfFile
print "Coord DCD file:", coordDcdFile

# start MDAnalysis with a universal
coord = MDAnalysis.Universe(inPsfFile,coordDcdFile)
dna = coord.selectAtoms("protein and name CA")
nAtoms = len(dna.atoms)
nSteps = len(coord.trajectory)

# print some log info
print "number of selected:", nAtoms
print "nSteps:", nSteps

#initialize average coordinates
avgCoord = numpy.zeros((nAtoms,3),dtype=float)
allCoord = numpy.zeros((nAtoms,3,nSteps),dtype=float)

# read in coordinates and add to average
for ts in coord.trajectory[0:nSteps]:
	# center the coordinates
	dna.translate(-dna.centerOfGeometry())
	# accumulate the average
	avgCoord += dna.positions	
	# populate the coordinate matrix
	allCoord[:,:,ts.frame-1] = dna.positions
#	if ts.frame%500==0:
#		print "Time step:", ts.frame

#complete the average
avgCoord /= float(nSteps)

# Reaverage to convergence
iteration = 0
residual = thresh + 10.0  # arbitrary assignment greater than thresh
newAvgCoord = numpy.zeros((nAtoms,3),dtype=float)
while residual > thresh and iteration < maxIter:
	newAvgCoord = numpy.zeros((nAtoms,3),dtype=float)
	# iterate over trajectory
	for ts in coord.trajectory[0:nSteps]:
#		dna.positions = avgCoord
#		dna.write("average_prealign.pdb")
#		dna.positions = allCoord[:,:,ts.frame-1]
#		dna.write("allcoord.pdb")
		# Compute rotation matrix to take time step coordinates and average coordinates
		R, d = rotation_matrix(allCoord[:,:,ts.frame-1],avgCoord)
		# apply rotation matrix
		allCoord[:,:,ts.frame-1] = numpy.dot(allCoord[:,:,ts.frame-1],R.T)
		# accumulate new average
		newAvgCoord += allCoord[:,:,ts.frame-1]

	# finish average
	newAvgCoord /= float(nSteps)
	# compute rmsd between new average and previous step average
	residual = MDAnalysis.analysis.align.rmsd(avgCoord,newAvgCoord)
	# update step info
	iteration += 1
	avgCoord = newAvgCoord
	# print log file
	print "Step:", iteration, "Residual", residual

print "Converged"

for i in range(0,nSteps):
	rmsd = MDAnalysis.analysis.align.rmsd(allCoord[:,:,i],avgCoord)
	print "%10d %10.5f" %(i+1, rmsd)
 
