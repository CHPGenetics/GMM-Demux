import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GMM_Demux",
    version="0.0.1",
    author="Hongyi Xin",
    author_email="gohongyi@gmail.com",
    description="A multiplet removal tool for processing cell hashing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CHPGenetics/GMM-demux",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
