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
        'Model page': 'https://csdms.colorado.edu/wiki/Model:OTTAR',
    },
)
