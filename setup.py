# usage: pip install --editable .
from setuptools import setup, find_packages
from openresim.__init__ import __version__
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

with open("requirements.txt") as f:
    requirements = f.read().splitlines()
    requirements = [r for r in requirements if r[:2] != "-e"]

setup(
    name="openresim",
    version=__version__,
    packages=find_packages(include=["openresim"]),
    author="ZAKARIYA ABUGRIN",
    author_email="zakariya.abugrin@gmail.com",
    description="Open Reservoir Simulation Library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="""
        Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
        International License
    """,
    install_requires=requirements,
    # setup_requires=['pytest-runner'],
    # tests_require=['pytest==4.4.1'],
    # test_suite='tests',
)
