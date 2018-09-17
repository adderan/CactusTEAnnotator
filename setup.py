from setuptools import setup, find_packages
import os
import subprocess


setup(
    name="CactusTEAnnotator",
    version="1.0",
    author="Alden Deran",
    packages=find_packages(where='.'),
    include_package_data=True,
    # We use the __file__ attribute so this package isn't zip_safe.
    zip_safe=False,

    install_requires=['networkx', 'multiset'],
    
    entry_points={
        'console_scripts': ['CactusTEAnnotator = CactusTEAnnotator.findRepeats:main', 'buildSubfamilies = CactusTEAnnotator.buildSubfamilies:main', 'treeBuilding = CactusTEAnnotator.treeBuilding:main', 'scoreGFF = CactusTEAnnotator.ari:main']})
