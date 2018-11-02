import shutil
import os

folders = [folder for folder in os.listdir('./') if os.path.isdir(folder)]
foldercont = []

for folder in folders:
	foldercont.extend(['./{0}/{1}'.format(folder, filev) for filev in os.listdir('./{0}'.format(folder)) if filev.endswith('.dat')])

os.mkdir('./copyfiles/')

for filev in foldercont:
	shutil.copy(filev, './copyfiles/{0}'.format(filev.split('/')[-1]))