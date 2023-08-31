# usage  (fixed): pip install .
# usage  (fixed ): python setup.py install
# usage (listen): pip install --editable .
# usage (listen): python setup.py develop
import re
from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


with open("requirements.txt") as f:
    requirements = f.read().splitlines()
    requirements = [r for r in requirements if r[:2] != "-e"]

with open("reservoirflow/__init__.py") as f:
    version = re.findall("__version__.*(\d.\d.\d).*", f.read())[0]

setup(
    name="reservoirflow",
    version=version,
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
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Reservoir Engineers",
        "Intended Audience :: Mathematicians",
        "Intended Audience :: Physicists",
        "Intended Audience :: Researchers",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    setup_requires=["pytest-runner"],
    tests_require=["unittest"],
    test_suite="tests",
)
