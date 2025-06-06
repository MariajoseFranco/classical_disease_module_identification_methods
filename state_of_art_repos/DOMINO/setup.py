from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name='domino-python',
    version="0.1.1",
    author="Hagai Levi",
    author_email="hagai.levi.007@gmail.com",
    description='DOMINO: Discovery of Modules In Networks using Omic',
    url='https://github.com/Shamir-Lab/DOMINO',
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    packages = find_packages(),
    package_data={'': ['*']},
    include_package_data=True,
    install_requires=[
        'networkx==2.4',
        'numpy==1.22.0',
        'scipy==1.10.0',
        'pandas==1.5.1',
        'pcst-fast==1.0.7',
        'statsmodels==0.11.0',
        'python-louvain==0.14'],
    entry_points = {
        "console_scripts": [
            "domino=src.runner:main_domino",
            "slicer=src.runner:main_slicer",
        ]
    }

)
