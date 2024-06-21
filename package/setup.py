import setuptools

setuptools.setup(

    name="pyctramer",   # The file name of the package after packaging 

    version="0.2.6",      # version number 

    author="Zengkui Liu, Dominikus Brian, and Xiang Sun",     # author name 

    author_email="zengkui.liu@nyu.edu",   # author email 

    description="Python package for Charge Transfer Rate from Atomistic Molecular dynamics, Electronic structure, and Rate theory",

    long_description="Python package for Charge Transfer Rate from Atomistic Molecular dynamics, Electronic structure, and Rate theory",

    long_description_content_type="text/markdown",  # Required dependencies

    install_requires=[],  

    url="https://www.github.com/ctramer/PyCTRAMER",  # Replace with your actual project URL 

    packages=setuptools.find_packages(),   # The python package directionry 

)
