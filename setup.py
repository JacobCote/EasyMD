from setuptools import setup, find_packages

setup(
    name="EasyMD",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "pyyaml",
        "mdtraj>=1.9.0",
    ],
    entry_points={
        "console_scripts": [
            "easymd = EasyMD.__main__:main"
        ]
    },
    include_package_data=True,
    zip_safe=False,
)