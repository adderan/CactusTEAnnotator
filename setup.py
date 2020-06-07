from setuptools import setup, find_packages
import os
import subprocess
from setuptools.command.install import install
from setuptools.command.develop import develop

class PostInstallCommand(install):
    def run(self):
        subprocess.run(["pip", "install", "cactus/submodules/sonLib"], check=True)
        install.run(self)

class PostDevelopCommand(develop):
    def run(self):
        subprocess.run(["pip", "install", "cactus/submodules/sonLib"], check=True)
        develop.run(self)

setup(
    name="CactusTEAnnotator",
    version="1.0",
    author="Alden Deran",
    packages=find_packages(where='.'),
    include_package_data=True,

    zip_safe=False,

    cmdclass = {
        'install': PostInstallCommand,
        'develop': PostDevelopCommand,
    },
   
    entry_points={
        'console_scripts': ['CactusTEAnnotator = CactusTEAnnotator.findRepeats:main', 'scoreGFF = CactusTEAnnotator.scoreGFF:main']}
)
