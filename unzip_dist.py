import tarfile
import sys
version = sys.argv[1]
tar = tarfile.open(f"dist/bioscrape-{version}.tar.gz")
for member in tar.getmembers():
    print(member)
tar.close()
