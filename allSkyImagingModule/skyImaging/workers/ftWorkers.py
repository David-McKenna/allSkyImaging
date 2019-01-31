"""Function that perform the processing maths behind the DFT/FFT methods
"""
import numpy as np

def dftWorker(idx, correlationMatrix, obsFreq, weight, conjWeight, skyView):
	"""DFT implementation for generating all sky images
	
	Args:
		idx (int): Multiprocessing reference ID
		correlationMatrix (np.ndarray): Array of antenna correlations
		obsFreq (float): Observing frequency (for debug messages)
		weight (np.ndarray): UVW plane weights
		conjWeight (np.ndarray): UVW conj. weights
		skyView (np.ndarray): Output array (for shape reference)

	Returns:
		idx, outputSkyImage: Multiprocessing ID, output results
	"""
	frameCount = correlationMatrix.shape[-1]
	for frame in np.arange(frameCount):
		correlationMatrixChan = correlationMatrix[..., frame] # (96,96)
		print("Processing Frame {0} of {1} in Chunk {2} at Frequency {3:.2F}E+06".format(frame + 1, frameCount, idx + 1, obsFreq / 1e6))

		tempProd = np.dot(conjWeight, correlationMatrixChan) # w.H * corr # (l,m, 1, 96)

		prodProd = np.multiply(tempProd.transpose((0,1,3,2)), weight)
		skyView[..., frame] = np.sum(prodProd, axis = (2,3))

	print("Chunk {0} completed work on {1} frames.".format(idx + 1, frameCount))
	return idx, skyView

def fftWorker(idx, correlationMatrix, uvIdx, convolvedWeight, skyView):
	"""FFT implementation for generating all sky images
	
	Args:
	    idx (int): Multiprocessing reference ID
	    correlationMatrix (np.ndarray): Array of antenna correlations
	    uvIdx (np.ndarray): UV plane indices
	    convolvedWeight (np.ndarray): UVW plane weights with FFT modifications
		skyView (np.ndarray): Output array (for shape reference)
	
	Returns:
		idx, outputSkyImage: Multiprocessing ID, output results

	"""
	print(uvIdx.shape, convolvedWeight.shape, correlationMatrix.shape)
	for idx, uvTuple in enumerate(zip(uvIdx[0, 0, :], uvIdx[0, 1, :])):
		skyView[uvTuple[0], uvTuple[0], :] += (convolvedWeight[:, idx, np.newaxis] * correlationMatrix[:, idx, :])[0, :]

	skyView[...] = np.fft.ifft2(np.fft.ifftshift(skyView, axes = (0,1)), axes = (0,1))
	skyView[...] = np.fft.ifftshift(skyView, axes = (0,1))

	return idx, skyView
