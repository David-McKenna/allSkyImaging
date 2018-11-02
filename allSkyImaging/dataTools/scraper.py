"""Summary
"""
import re
import urllib2
import urllib

from bs4 import BeautifulSoup

def scrape(url, outputFolder = './', dataFiles = True, subbandScan = False):
	"""Summary
	
	Args:
	    url (TYPE): Description
	    outputFolder (str, optional): Description
	    dataFiles (bool, optional): Description
	    subbandScan (bool, optional): Description
	
	Returns:
	    TYPE: Description
	"""
	sbReg = re.compile(r'sb\d*')
	datReg = re.compile(r'(?=\.dat).*')
	
	soup = [BeautifulSoup(urllib2.urlopen(url))]

	if subbandScan:
		soup = soup[0]
		sbLinks = soup.find_all(href=sbReg)
		soup = [BeautifulSoup(urllib2.urlopen(url + link['href'])) for link in sbLink]
		links = [url + link['href'] for link in sbLink]

	else:
		links = [url]


	os.mkdir(outputFolder)
	scrapedFiles = []
	for idx, soupVar in enumerate(soup):
		dataLinks = soupVar.find_all(href=dataReg)

		for link in dataLinks:
			print(links[idx] + link['href'])

			ouputFile = outputFolder + link['href']
			urllib.urlretrieve(links[idx] + link['href'], outputFile)
			scrapedFiles.append(outputFile)

	print('Scrape Completed')

	return scrapedFiles