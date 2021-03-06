import sys
from setuptools import setup
import versioneer

if sys.version_info < (3, 0):
	sys.stdout.write("At least Python 3.0 is required.\n")
	sys.exit(1)

with open('README.md') as f:
	long_description = f.read()


setup(
	name = 'secapr',
	version = versioneer.get_version(),
	cmdclass = versioneer.get_cmdclass(),
	author = 'Tobias Andermann',
	author_email = 'tobias.andermann@bioenv.gu.se',
	url = 'https://github.com/AntonelliLab/seqcap_processor',
	description = 'Process sequence-capture fastq files into alignments for phylogenetic analyses',
	long_description = long_description,
	license = 'MIT',
	entry_points = {'console_scripts': ['secapr = secapr.__main__:main']},
	packages = ['secapr'],
    install_requires = [
	    # No dependencies listed here since we need to rely on conda anyway
    ],
	classifiers = [
		"Development Status :: 4 - Beta",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 2.7",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
