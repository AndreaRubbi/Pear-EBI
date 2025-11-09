from pathlib import Path
import os
import subprocess
import sys

from setuptools import Extension, find_packages, setup
from setuptools.command.build_py import build_py as _build_py
import shutil


def build_native_tools(package_root: str):
    """Attempt to build native tools (HashRF, tqDist) in the package during
    the build step. This runs `make` in HashRF and runs CMake/make in tqDist if
    available. If these tools are not available the commands will raise and the
    build will fail (which is appropriate when a platform-native binary must be
    compiled).

    Args:
        package_root: path to package directory containing `calculate_distances`.
    """
    cd = os.getcwd()
    try:
        calc_dir = os.path.join(package_root, "calculate_distances")

        # Build HashRF (has a Makefile)
        hashrf_dir = os.path.join(calc_dir, "HashRF")
        if os.path.exists(hashrf_dir):
            print(f"Building native HashRF in {hashrf_dir}")
            # run make; allow parallelism
            try:
                subprocess.run(["make", "-C", hashrf_dir], check=True)
            except subprocess.CalledProcessError as e:
                print(
                    "Warning: Failed to build HashRF. Ensure build tools (make, gcc/clang) are installed if you need the HashRF executable."
                )
                # continue without raising; wheel can still be built and include sources

        # Build tqDist (uses CMake)
        tq_dir = os.path.join(calc_dir, "tqDist")
        if os.path.exists(tq_dir):
            print(f"Building native tqDist in {tq_dir}")
            # Prefer using the package helper if available for a clearer
            # cross-platform build flow and better error messages.
            try:
                # Import from source tree; this is safe during sdist/wheel builds
                import pear_ebi._install_helpers as _ih

                ok, message = _ih.build_tqdist(package_root)
                if ok:
                    print("tqDist build succeeded:", message)
                else:
                    print("tqDist build helper reported failure:", message)
                    # fall back to the legacy steps below
                    raise RuntimeError(message)
            except Exception:
                # Legacy fallback: try to configure/build with cmake or make.
                build_subdir = os.path.join(tq_dir, "build")
                os.makedirs(build_subdir, exist_ok=True)
                # prefer cmake if available
                cmake = shutil.which("cmake")
                if cmake:
                    try:
                        subprocess.run([cmake, ".."], cwd=build_subdir, check=True)
                        subprocess.run(["cmake", "--build", "."], cwd=build_subdir, check=True)
                    except subprocess.CalledProcessError:
                        print("Warning: CMake build for tqDist failed. Will try plain make as fallback.")
                        try:
                            subprocess.run(["make", "-C", tq_dir], check=True)
                        except subprocess.CalledProcessError:
                            print(
                                "Warning: Failed to build tqDist. Install CMake or make and a C++ compiler if you need the tqDist executables."
                            )
                            # continue without raising
                else:
                    # fallback to running make in tqDist
                    try:
                        subprocess.run(["make", "-C", tq_dir], check=True)
                    except subprocess.CalledProcessError:
                        print(
                            "Warning: Failed to build tqDist (no cmake available). Install cmake or make + compiler if you need the tqDist executables."
                        )
                        # continue without raising

    finally:
        os.chdir(cd)


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
    version="1.0.1.5",
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
