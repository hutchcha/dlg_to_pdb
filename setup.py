from setuptools import setup, find_packages

setup(
    name="dlg_to_pdb",  # The package name
    version="0.1.0",  # Initial version
    author="Chase Hutchins",  # Replace with your name
    author_email="chasehutchins@comcast.net",  # Replace with your email
    description="A Python tool to convert AutoDock DLG files into PDB files",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",  # Format for README
    url="https://github.com/hutchcha/dlg_to_pdb",  # Replace with your GitHub repo URL
    packages=find_packages(),  # Automatically find package directories
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",  # Specify Python version
    entry_points={
        "console_scripts": [
            "dlg_to_pdb=dlg_to_pdb.dlg_to_pdb:main",  # Command-line entry point
        ],
    },
    install_requires=[
        # List any dependencies here, for example:
        # "numpy>=1.21.0"
    ],
    include_package_data=True,  # Include other files (e.g., README, LICENSE)
)
