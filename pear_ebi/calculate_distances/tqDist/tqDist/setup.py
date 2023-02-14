from distutils.core import setup
from platform import system

if system() == "Windows":
	ext = ".dll"
else:
	ext = ".so"

setup(name = "tqDist",
      version = "1.0.0",
      description = "Package for computing the triplet or quartet distances between rooted or unrooted trees (binary or arbitrary degree).",
      author = "Andreas Sand",
      author_email = "asand@birc.au.dk",
      url = "http://birc.au.dk/software/tqDist",
      packages = ["pyTQDist"],
      package_data = {"": ["libpyQuartetDist%s" % ext, "libpyTripletDist%s" % ext]},
      license = "LGPL" 
    )
