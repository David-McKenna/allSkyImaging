"""data.lofar.ie scraper
"""
import re
import os

# Why do I need to use both to get this to function properly...
import urllib2
import urllib

from bs4 import BeautifulSoup

def scrape(url, outputFolder = './', subbandScan = False):
	"""Scrape data.lofar.ie for when you're stuck on Trinity's god awful network but need data.
	
	Args:
	    url (str): Folder of observation
	    outputFolder (str, optional): Where to dump the data
	    subbandScan (bool, optional): If we ran a subband scan, scrpae the subbfolders.
	
	Returns:
	    list: List of scraped folders
	"""
	sbReg = re.compile(r'sb\d*')
	dataReg = re.compile(r'(?=\.dat).*')
	
	soup = [BeautifulSoup(urllib2.urlopen(url))]

	if subbandScan:
		soup = soup[0]
		sbLinks = soup.find_all(href=sbReg)
		soup = [BeautifulSoup(urllib2.urlopen(url + link['href'])) for link in sbLinks]
		links = [url + link['href'] for link in sbLinks]

	else:
		links = [url]


	os.makedirs(outputFolder)
	scrapedFiles = []
	for idx, soupVar in enumerate(soup):
		dataLinks = soupVar.find_all(href=dataReg)

		for link in dataLinks:
			print(links[idx] + link['href'])

			outputFile = outputFolder + link['href']
			urllib.urlretrieve(links[idx] + link['href'], outputFile)
			scrapedFiles.append(outputFile)

	print('Scrape Completed')

	return scrapedFiles
	