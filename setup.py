import setuptools

# Read README for long description
with open("README.md", "r", encoding="utf-8") as readme_file:
    LONG_DESCRIPTION = readme_file.read()

setuptools.setup(
    name="MLDockKit",
    version="0.1.1",
    include_package_data=True,
    description="Python package that calculates Lipinski descriptors, predicts pIC50, and performs docking",
    url="https://github.com/clabe-wekesa/MLDockKit",
    license="MIT",
    packages=setuptools.find_packages(),
    install_requires=[
        "rdkit",
        "pandas",
        "padelpy",
        "joblib",
        "scikit-learn",  # Removed incorrect version constraint
        "scipy",
        "gemmi",
        "numpy",
    ],
    python_requires=">=3.9",
    author="Edwin Mwakio, Clabe Wekesa, Patrick Okoth",
    author_email="simiyu86wekesa@gmail.com",
    package_data={
        # Include your non-Python files here
        "MLDockKit": [
            "padel_model.joblib",
            "7te7_H_no_HOHandI0V.pdb",
            "7te7_original.pdb",
            "7te7_prepared.pdbqt",
        ],
    },
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    classifiers=[
        "Natural Language :: English",
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
