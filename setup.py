from deeprho import __version__
from setuptools import setup, find_packages

CLASSIFIERS = [
    'Development Status :: 1',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: Apache License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Topic :: Population Genetics',
    'Operating System :: POSIX :: Linux',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']


setup(
    name="deeprho",
    version=__version__,
    author="Haotian Zhang, Yufeng Wu",
    author_email="haotianzh@uconn.edu, yufeng.wu@uconn.edu",
    long_description="deeprho is a software used for estimating recombination rate given population genetic data",
    url="https://github.com/haotianzh/deeprho_v2",
    packages=find_packages(),
    package_data={
        "deeprho.popgen": ["libs/*.jar"]
    },
    entry_points={
        "console_scripts": ["deeprho = deeprho.entry:main"]
    },
    python_requires=">=3.7",
    # setup_requires=['setuptools==57.5.0'],
    install_requires= ["h5py>=2.10.0",
                       "JPype1==1.3.0",
                       "matplotlib>=3.5.1",
                       "pandas>=1.0.5",
                       "numpy>=1.19.2",
                       "pptree==3.1",
                       "scikit-learn>=0.23.1",
                       "scipy>=1.6.2",
                       "tskit>=0.3.5",
                       "coloredlogs"
                       ],
    extras_require={"simulate": ["msprime>=1.0.0"], "estimate":["pyvcf", "tensorflow"]},
    keywords= "population genetics, recombination, algorithm",
    classifiers=CLASSIFIERS
)



