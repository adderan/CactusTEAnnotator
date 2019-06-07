from setuptools import setup, find_packages
import os
import subprocess


setup(
    name="CactusTEAnnotator",
    version="1.0",
    author="Alden Deran",
    packages=find_packages(where='.'),
    include_package_data=True,

    zip_safe=False,

    install_requires=['actualSonLib'],
    
    entry_points={
        'console_scripts': ['CactusTEAnnotator = CactusTEAnnotator.findRepeats:main', 'scoreGFF = CactusTEAnnotator.scoreGFF:main']})
