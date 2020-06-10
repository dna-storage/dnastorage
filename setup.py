import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dnastorage", 
    version="0.9.0",
    author="James M. Tuck",
    author_email="james.m.tuck@gmail.com",
    description="dnastorage provides basic support for encoding and decoding files into and out of DNA with several convenient interfaces.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dna-storage/dnastorage",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
