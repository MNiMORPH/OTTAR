import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="OTTAR",
    version="0.2.0",
    author="Andrew D. Wickert",
    author_email="awickert@umn.edu",
    description="Ode To Transient (Ancho de los) Rivers: Transient evolution of river-channel width",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MNiMORPH/ottar",
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'pandas',
        ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Intended Audience :: Science/Research",
    ],
    keywords='fluvial geomorphology sediment transport landscape evolution',
    project_urls={
        'Source and README': 'https://github.com/MNiMORPH/OTTAR',
        'CSDMS repository': 'https://csdms.colorado.edu/wiki/Model:OTTAR',
        'Zenodo': 'https://zenodo.org/record/5781792/export/dcat#.Yme_nIyxXCK',
        'DOI': 'https://dx.doi.org/10.5281/zenodo.5124965',
    },
)
