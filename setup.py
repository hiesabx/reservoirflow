# pip install --editable .

from setuptools import setup, find_packages

setup(
    name = 'openresim',
    packages = find_packages(include=['openresim']),
    version='0.0.1',
    description='The Pythonic Petroleum Reservoir Simulator',
    author='Zakariya AbuGrin',
    license='MIT',
    # install_requires=[],
    # setup_requires=['pytest-runner'],
    # tests_require=['pytest==4.4.1'],
    # test_suite='tests',
)