import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SCYN",
    version="1.0.5",
    author="Xikang Feng",
    author_email="xikanfeng2@gmail.com",
    description="SCYN: Single cell CNV profiling method using dynamic programming",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xikanfeng2/SCYN",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
