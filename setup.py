from setuptools import setup, find_packages

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='trilogy',
    version='0.0.1',
    license='GPLv3',
    author="Dan Coe",
    author_email='dcoe@stsci.edu',
    maintainer='Renan Alves de Oliveira',
    maintainer_email='fisica.renan@gmail.com',
    description='Python script that converts astronomical FITS images in color/grayscale images.',
    url='https://github.com/oliveirara/trilogy',
    keywords='fits files,astronomical images',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['trilogy'],
    scripts=['trilogy/trilogy-cl'],
    install_requires=[
        'Pillow>=8','astropy>=4','numpy>=1.16','scipy>=1.7',
    ],
)

