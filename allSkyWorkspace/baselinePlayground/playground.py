import numpy as np
import matplotlib as mpl

#mpl.use('Agg')

import matplotlib.pyplot as plt
import h5py
import datetime
import logging
import itertools
import collections
from numba import jit


def initArray():
	antLoc = np.load('antLoc.npy')
	hbaDelta = np.load('hbadeltas.npy')

	antPos = antLoc[..., 0] # Only consider the xpolarisation... theyre all the same anyway.

	return antPos, hbaDelta

#@jit(nopython=True)
def costFunc(arrayPos):
	baselines = arrayPos[..., np.newaxis] - arrayPos.T[np.newaxis, ...]
	#deltaSum = -1. * np.sum(np.abs(baselines)) / (5. + 0.15)
	deltaBaselines = -1. * np.sum(np.abs(baselines[..., np.newaxis, np.newaxis] - baselines[np.newaxis, np.newaxis, ...])) / (5.15)

	return deltaBaselines

#@jit(nopython=True)
def costFuncDFT(beamIm):
	beamShape = beamIm.shape
	initMax = np.max(beamIm)
	beamIm[beamIm > initMax / (2. * np.e) ] = 0.

	postMax = np.max(beamIm)
	beamIm[beamIm > postMax / (2. * np.e) ] = 0.

	finalMax = np.max(beamIm)
	finalAvg = np.mean(beamIm)
	return postMax + finalMax + finalAvg - initMax

global tempConst
tempConst = collections.deque([1e16], maxlen = 128)
def probFlip(deltaSumOrig, deltaSumNew, currLen):
	deltaDiff = (deltaSumNew - deltaSumOrig)
	tempAvg = 1.
	expTerm = np.exp(-1. * deltaDiff / tempAvg)

	tempConst.append(deltaDiff)

	print(expTerm)
	return expTerm

def heartbeat(seed = 0):
	np.random.seed(seed)
	antPos, hbaDelta = initArray()

	#processes = mp.cpu_count() - 1
	#mpPool = mp.Pool(processes = processes)

	i = 0
	#__ = [mpPool.apply_async(initThread, args = ([i, seed + i, antPos, hbaDelta])) for i in range(processes)]
	initThread(i, seed + i, antPos, hbaDelta)

#@jit(nopython=True)
def initThread(i, seed, antPos, hbaDelta):
	np.random.seed(seed)
	with h5py.File('./{0}-{1}.h5'.format(datetime.datetime.utcnow(), i), 'a') as checkpointRef:
		stepCount = 0
		arraySize = 96
		antChoice = 16
		activationChoice = np.random.randint(0, 16, (arraySize,), dtype = int)

		arrayLayout = antPos + hbaDelta[activationChoice]

		stepCheckpoint = checkpointRef.create_dataset(str(stepCount), (arraySize,), compression = 'lzf', dtype = np.int8)
		stepCheckpoint[...] = activationChoice

		prevCost = costFunc(arrayLayout)
		#globalMin = 0
		prevCost = 1e16
		globalMin = 1e16
		stepCheckpoint.attrs.create('cost', prevCost)

		stepsPerLoop = 256
		changesPerStep = 16
		while True:
			activationSwap = np.random.randint(0, 16, (stepsPerLoop, changesPerStep))
			antSelect = np.random.randint(0, 96, (stepsPerLoop, changesPerStep))
			refRand = np.random.rand(1)

			for actIdx, antIdx in itertools.izip(activationSwap, antSelect):
				print(actIdx, antIdx, stepCount, prevCost)
				stepCount += 1
				cacheAnt, cacheAct = arrayLayout[antIdx], activationChoice[antIdx]

				arrayLayout[antIdx] = antPos[antIdx] + hbaDelta[actIdx]
				activationChoice[antIdx] = actIdx

				#newCost = costFunc(arrayLayout)

				beamIm = dftImage(arrayLayout)
				newCost = costFuncDFT(beamIm)

				print(newCost, prevCost)

				if newCost < prevCost:
					print('Changed 0')
					prevCost = newCost

					if prevCost < globalMin:
						globalMin = prevCost
						stepCheckpoint = checkpointRef.create_dataset(str(stepCount), (arraySize,), compression = 'lzf', dtype = np.int8)
						stepCheckpoint[...] = activationChoice
						stepCheckpoint.attrs.create('cost', globalMin)
					continue
				else:
					if probFlip(prevCost, newCost, changesPerStep) > np.random.rand(1):
						print('Changed 1')
						prevCost = newCost
						continue
					else:
						print('Skipped')
						arrayLayout[antIdx] = cacheAnt
						activationChoice[antIdx] = cacheAct

			stepCheckpoint = checkpointRef.create_dataset(str(stepCount), (arraySize,), compression = 'lzf', dtype = np.int8)
			stepCheckpoint[...] = activationChoice

			prevCost = costFunc(arrayLayout)
			stepCheckpoint.attrs.create('cost', prevCost)
			changesPerStep -= 4
			changesPerStep = max(changesPerStep, 1)

			if stepCount > 100000:
				return

#@jit(nopython=True, parallel=True)
def dftImage(uvw,px=[100, 100],res=100 / np.pi,mask=False):
	"""return a DFT image"""
	uvw = np.array(getBaselines(uvw))
	nants=uvw.shape[0]
	
	im=np.zeros((px[0],px[1]),dtype=complex)

	mid_k=int(px[0]/2.)
	mid_l=int(px[1]/2.)

	u=uvw[:,0]
	v=uvw[:,1]

	u/=mid_k
	v/=mid_l

	# Speedup with numpy arrays / numba port in future
	im = jitDFT(im, px, u, v, mid_k, mid_l)

	return im.real

@jit(nopython=True, parallel=True)
def jitDFT(im, px, u, v, mid_k, mid_l):
	for k in range(px[0]):
		for l in range(px[1]):
			im[k,l] = np.sum(np.exp(-2.*np.pi*1j*((u*(k-mid_k)) + (v*(l-mid_l)))))
	return im

#@jit(nopython=True,parallel=True)
def getBaselines(antLoc):
	baselines = []
	for idx in range(len(antLoc)):
		ant = antLoc[idx]
		for idx2 in range(len(antLoc[idx:])):
			diff = ant - antLoc[idx + idx2]
			baselines.append(diff)
			baselines.append(-1. * diff)

	return baselines