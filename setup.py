import setuptools

# Read README for long description
with open('README.md', 'r', encoding='utf-8') as readme_file:
	LONG_DESCRIPTION = readme_file.read()

setuptools.setup(
    name="MLDockKit",
    version="0.0.1",
    include_package_data=True,
    description='Python package that calculates Lipinsky descriptors, predicts pIC50 and performs docking',
    url = "https://github.com/clabe-wekesa/MLDockKit",
    license = "MIT",
    packages=setuptools.find_packages(),
    install_requires=[
        # List your Python dependencies here
        'rdkit', 'pandas', 'padelpy', 'joblib', 'meeko', 'scikit-learn', 'scipy', 'numpy'
    ],
    author='Edwin Mwakio, Clabe Wekesa, Patrick okoth',
    author_email='simiyu86wekesa@gmail.com',  
    package_data={
        # Include your non-Python files here
        'MLDockKit': ['padel_model.joblib', '5gs4.pdbqt'],
    },
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    classifiers = [
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Operating System :: OS Independent"
],
)

