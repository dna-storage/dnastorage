import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dnastorage", 
    version="0.9.1",
    author="James M. Tuck",
    author_email="james.m.tuck@gmail.com",
    description="dnastorage provides basic support for encoding and decoding files into and out of DNA with several convenient interfaces.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dna-storage/dnastorage",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
    ],
    python_requires='>=3.6',
)
