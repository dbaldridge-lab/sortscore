from setuptools import setup, find_packages

# Read requirements from requirements.txt
def read_requirements():
    with open('requirements.txt', 'r') as f:
        lines = f.readlines()
    # Filter out comments and empty lines
    requirements = []
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#'):
            requirements.append(line)
    return requirements

setup(
    name='sortscore',
    version='0.1.0',
    author='Caitlyn Chitwood',
    author_email='c.chitwood@wustl.edu',
    description='A modular Python package for Sort-seq variant analysis',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/dbaldridge-lab/sortscore',
    packages=find_packages(include=["sortscore", "sortscore.*"]),
    install_requires=read_requirements(),
    python_requires='>=3.10',
    entry_points={
        'console_scripts': [
            'sortscore=sortscore.run_analysis:main',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    keywords='bioinformatics sequencing variant-analysis sort-seq',
    include_package_data=True,
    zip_safe=False,
)