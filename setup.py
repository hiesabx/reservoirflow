# usage  (fixed): pip install .
# usage  (fixed ): python setup.py install
# usage (listen): pip install --editable .
# usage (listen): python setup.py develop

from setuptools import setup, find_packages
from reservoirflow.__init__ import __version__
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

with open("requirements.txt") as f:
    requirements = f.read().splitlines()
    requirements = [r for r in requirements if r[:2] != "-e"]

setup(
    name="reservoirflow",
    version=__version__,
    author="Zakariya Abugrin",
    author_email="zakariya.abugrin@gmail.com",
    maintainer="Zakariya Abugrin",
    maintainer_email="zakariya.abugrin@gmail.com",
    description="a Petroleum Reservoir Simulation Library in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zakgrin/reservoirflow",
    download_url="https://github.com/zakgrin/reservoirflow.git",
    license="""
        Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
        International License
    """,
    keywords=["Petroleum", "Reservoir", "Simulation", "Scientific Computing"],
    project_urls={
        "numpy": "https://numpy.org/",
        "scipy": "https://scipy.org/",
        "sympy": "https://sympy.org/",
    },
    python_requires=">=3.7",
    packages=find_packages(include=["reservoirflow"]),
    include_package_data=True,
    install_requires=requirements,
    classifiers=[
        "Intended Audience :: Reservoir Engineers",
        "Intended Audience :: Researchers",
        "Intended Audience :: Academics",
        "Intended Audience :: Students",
        "Programming Language :: Python :: 3",
    ],
    setup_requires=["pytest-runner"],
    tests_require=["unittest"],
    test_suite="tests",
)
