from pathlib import Path

from setuptools import find_packages, setup

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="pear_ebi",
    version="0.1.1.10",
    license="MIT License",
    description="Embeds phylogenetic tree distances and produce representations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Andrea Rubbi",
    author_email="andrea.rubbi.98@gmail.com",
    url="https://github.com/AndreaRubbi/TreeEmbedding",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "scikit-learn",
        "plotly",
        "rich",
    ],
)
