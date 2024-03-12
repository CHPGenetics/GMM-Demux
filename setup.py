import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GMM_Demux",
    version="0.2.2.3",
    author="Hongyi Xin",
    author_email="gohongyi@gmail.com",
    description="A multiplet removal tool for processing cell hashing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CHPGenetics/GMM-demux",
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'GMM-demux = GMM_Demux.GMM_Demux:main'
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "pandas>=1.4.3",
        "numpy>=1.22.4",
        "scipy>=1.12.0",
        "tabulate",
        "argparse",
        "statistics",
        "BitVector",
        "scikit-learn"
    ]
)
