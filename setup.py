from setuptools import setup, find_packages

setup(
    name='crypsplice',
    version='2.0.0',
    description='CrypSplice 2.0: Detect and quantify cryptic splicing events from RNA-seq data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Michelle Dias',
    author_email='michelledias10@gmail.com',
    url='https://github.com/michelle-dias/CrypSplice2.0',  # <-- update this
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    install_requires=open('requirements.txt').read().splitlines(),
    entry_points={
        'console_scripts': [
            'crypsplice=crypsplice.main:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
