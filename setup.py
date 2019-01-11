import os
import sys

if len(sys.argv) == 1:
	print('Installing Bioscrape')
	os.system("python setup_bioscrape.py install")
	print("Installing Bioscrape Lineages")
	os.system("python setup_lineages.py install")
elif len(sys.argv) == 2:
	if sys.argv[1] == "install":
		print('Installing Bioscrape')
		os.system("python setup_bioscrape.py install")
		print("Installing Bioscrape Lineages")
		os.system("python setup_lineages.py install")
	elif sys.argv[1] == "bioscrape":
		print('Installing Bioscrape')
		os.system("python setup_bioscrape.py install")
	elif sys.argv[1] == "lineages":
		print('Installing lineages')
		os.system("python setup_lineages.py install")
	else:
		raise ValueError("setup.py accepts the arguments: install, bioscrape, and lineages. By default, will install both bioscrape and lineages unless just one is specified.")
elif len(sys.argv) == 3:
	if sys.argv[1] != "install":
		raise ValueError("setup.py accepts the arguments: install, bioscrape, and lineages. By default, will install both bioscrape and lineages unless just one is specified.")
	else:
		if sys.argv[2] == "bioscrape":
			print('Installing Bioscrape')
			os.system("python setup_bioscrape.py install")
		elif sys.argv[2] == "lineages":
			print('Installing lineages')
			os.system("python setup_lineages.py install")
		else:
			raise ValueError("setup.py accepts the arguments: install, bioscrape, and lineages. By default, will install both bioscrape and lineages unless just one is specified.")