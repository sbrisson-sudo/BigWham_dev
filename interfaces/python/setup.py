from setuptools import setup, find_packages

setup(
    name='bigwham4py',
    version='0.1.0',
    packages=find_packages(),
    package_data={
        'bigwham4py': ['*.py', '*.so']
    },
    include_package_data=True,
)