from setuptools import setup, find_packages

with open("readme.md", "r") as fh:
    long_description = fh.read()

setup(
    name="mutperiodpy",
    version="0.7.0",
    description='Analyze mutational periodicity about nucleosomes',
    long_description_content_type="text/markdown",
    url='https://github.com/bmorledge-hampton19/Nucleosome-Mutation-Analysis',
    author='Ben Morledge-Hampton',
    author_email='b.morledge-hampton@wsu.edu',
    license='MIT',
    python_requires='>=3.7',
    packages=find_packages(),
    package_data={"mutperiodpy": ["run_mutperiodR/*.r", "run_mutperiodR/*.R", "Tkinter_scripts/test_tube.png",
                  "input_parsing/*.tsv"]},
    entry_points=dict(
        console_scripts=['mutperiod=mutperiodpy.Main:main']
    ),
    install_requires=["benbiohelpers"]
    
)