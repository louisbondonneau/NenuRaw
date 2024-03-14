from setuptools import setup, find_packages

setup(
    name='NenuRaw',
    version='0.01',
    python_requires='>=3.8',
    description='this package can handle RAW observation files of the radio telescopes of Nançay',
    author='Louis Bondonneau',
    author_email='louis.bondonneau@obs-nancay.fr',
    # packages=find_packages(),
    packages=find_packages(),
    package_data={"demo": ["*.py"],
                 "beta_dir": ["*.py"]},
    install_requires=[
        'numpy',
        'astropy',
		'matplotlib',
    ],
    scripts=["main.py"],
)
