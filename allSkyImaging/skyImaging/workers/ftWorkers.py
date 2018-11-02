"""Summary
"""

def dftWorker(idx, correlationMatrix, obsFreq, conjWeight, skyView):
	"""Summary
	
	Returns:
	    TYPE: Description
	
	Args:
	    idx (TYPE): Description
	    correlationMatrix (TYPE): Description
	    obsFreq (TYPE): Description
	    conjWeight (TYPE): Description
	    skyView (TYPE): Description
	"""
	frameCount = correlationMatrix.shape[2]
	for frame in np.arange(frameCount):
		correlationMatrixChan = correlationMatrix[..., frame] # (96,96)
		print("Processing Frame {0} of {1} at Frequency {2:.2F}E+06".format(frame + 1, frameCount, obsFreq / 1e6))

		tempProd = np.dot(conjWeight, correlationMatrixChan) # w.H * corr # (l,m, 1, 96)

		prodProd = np.multiply(tempProd.transpose((0,1,3,2)), weight).real
		skyView[..., frame] = np.sum(prodProd, axis = (2,3))

	return idx, skyView

def fftWorker(idx, correlationMatrix, uv, kernel, skyView):
	"""Summary
	
	Args:
	    idx (TYPE): Description
	    correlationMatrix (TYPE): Description
	    uv (TYPE): Description
	    kernel (TYPE): Description
	    skyView (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	if 'gaus' in kernel:
		gaussCoeff = float(kernel.split('/')[-1])

		pixelBoundary = int((skyView.shape[0] / 1.05) * 0.05)
		sampleArea = min(max(int(gaussCoeff * 10), 1), pixelBoundary) # Samples an area effectively ~5x the length/width of the Gaussian width

		print(sampleArea, gaussCoeff, pixelBoundary)
		guassMatrix = np.square(np.array([[gaussCoeff, 0.], [0., gaussCoeff]]))
		convKernel = lambda sampledPoints: (1. / np.sqrt(8. * np.pi) / gaussCoeff ) * np.sum(np.exp(-1. * np.dot(guassMatrix, np.square(sampledPoints))), axis = (0))

	else:
		convKernel = lambda sampledPoints: np.full(sampledPoints.shape[2:], fill_value = 1. / len(sampledPoints.shape[2:]))

	if 'rect' in kernel:
			sampleArea = 1

	uv = uv.astype(float)
	uv *= 0.967
	sampleCache = np.mgrid[-sampleArea + 1:sampleArea, -sampleArea + 1:sampleArea].reshape(2, -1).astype(float).T[..., np.newaxis] # (nSamples, 2, 1)


	offsets = (1. - (uv % 1.))[np.newaxis] # (1, 2, 96 * 96)

	sampleOffsets = sampleCache + offsets # (nSamples, 2, 96 * 96)
	uvIdx = (sampleOffsets + uv).astype(int)


	convolvedWeight = convKernel(sampleOffsets) # (nSamples, 96 * 96)
	print(sampleCache.shape, offsets.shape, uvIdx.shape, convolvedWeight.shape)
	plt.imshow(convolvedWeight[:, 0].reshape(sampleArea * 2 + 1, sampleArea * 2 + 1))
	plt.show()
	print(convolvedWeight)
	skyView[uvIdx[:, 0], uvIdx[:, 1], :] += convolvedWeight[..., np.newaxis] * correlationMatrix

	plt.imshow(skyView[..., 0])
	plt.show()
	skyView[...] = np.fft.ifft2(np.fft.ifftshift(skyView))
	skyView[...] = np.fft.ifftshift(np.abs(skyView.real)) # Abs of FFT == DFT?

	return idx, skyView
