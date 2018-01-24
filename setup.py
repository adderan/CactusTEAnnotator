from setuptools import setup, find_packages
import os
import subprocess


setup(
    name="repeatAnnotator",
    version="1.0",
    author="Alden Deran",
    package_dir = {'': 'src'},
    packages=find_packages(where='src'),
    include_package_data=True,
    # We use the __file__ attribute so this package isn't zip_safe.
    zip_safe=False,

    install_requires=[],
    
    entry_points={
        'console_scripts': ['repeatAnnotator = RepeatAnnotator.repeatAnnotator:main']},)
