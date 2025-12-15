from pathlib import Path
import os
import subprocess
import sys

from setuptools import Extension, find_packages, setup
from setuptools.command.build_py import build_py as _build_py
import shutil


def build_native_tools(package_root: str):
    """
    Ensure all binaries in linux_bin/HashRF, linux_bin/tqDist, mac_bin/HashRF, mac_bin/tqDist are executable.
    This is the only install-time action; no moving or copying is performed.
    """
    for folder in [
        os.path.join(package_root, "calculate_distances", "linux_bin", "HashRF"),
        os.path.join(package_root, "calculate_distances", "linux_bin", "tqDist"),
        os.path.join(package_root, "calculate_distances", "mac_bin", "HashRF"),
        os.path.join(package_root, "calculate_distances", "mac_bin", "tqDist"),
    ]:
        if os.path.isdir(folder):
            for root, dirs, files in os.walk(folder):
                for fname in files:
                    fpath = os.path.join(root, fname)
                    try:
                        os.chmod(fpath, 0o755)
                    except Exception as e:
                        print(f"Warning: could not set executable permission for {fpath}: {e}")



class build_py(_build_py):
    def run(self):
        # Ensure native tools are compiled before Python build copies package data
        package_root = os.path.join(os.path.dirname(__file__), "pear_ebi")
        try:
            build_native_tools(package_root)
        except Exception as e:
            print("Native build failed:", e)
            # Propagate exception so build/wheel creation fails — this is
            # preferable to shipping a wheel without usable native binaries.
            raise
        super().run()

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name="pear_ebi",
    version="1.0.1.6",
    license="MIT License",
    description="Embeds phylogenetic tree distances and produce representations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Andrea Rubbi",
    author_email="andrea.rubbi.98@gmail.com",
    url="https://github.com/AndreaRubbi/TreeEmbedding",
    packages=find_packages(),
    include_package_data=True,
    # Ensure non-Python package data (native binaries and helper files) are
    # included in built wheels/sdists. This makes files under
    # pear_ebi/calculate_distances/HashRF and pear_ebi/calculate_distances/tqDist
    # available at runtime after pip install.
    package_data={
        "pear_ebi": [
            "calculate_distances/HashRF/*",
            "calculate_distances/HashRF/*/*",
            "calculate_distances/tqDist/*",
            "calculate_distances/tqDist/*/*",
            "calculate_distances/linux_bin/*",
            "calculate_distances/linux_bin/*/*",
            "calculate_distances/mac_bin/*",
            "calculate_distances/mac_bin/*/*",
        ]
    },
    cmdclass={
        "build_py": build_py,
    },
    zip_safe=False,
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 3 - Alpha",
        # "Intended Audience :: Bioinformaticians",
        # "Topic :: Phylogenetics",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    install_requires=[
        "numpy<2.0.0",
        "pandas<=2.3.3",
        "matplotlib>=3.4",
        "scikit-learn<=1.6.1",
        "rich<=14.2.0",
        "pyDRMetrics==0.0.7",
        "tqdm<=4.67.1",
        "toml==0.10.2",
        "kaleido==1.2.0",
        "ipykernel==6.17.1",
        "ipython==8.6.0",
        "ipywidgets==7.7.2",
        "jupyter==1.0.0",
        "notebook==6.5.6",
        "jupyterlab==3.5.0",
        "nbconvert==6.4.5",
        "pandoc==2.4",
        "plotly==5.11.0",
        "scipy<=1.13.1",
        "Wand==0.6.13",
    ],
    entry_points={"console_scripts": ["pear_ebi = pear_ebi.__main__:main"]},
)
