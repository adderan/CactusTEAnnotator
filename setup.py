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

    install_requires=['actualSonLib', 'networkx', 'multiset'],
    
    entry_points={
        'console_scripts': ['CactusTEAnnotator = CactusTEAnnotator.findRepeats:main', 'scoreGFF = CactusTEAnnotator.ari:main']})
