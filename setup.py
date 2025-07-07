from setuptools import setup, find_packages

setup(
    name='sortscore',
    version='0.1.0',
    author='Caitlyn Chitwood',
    author_email='c.chitwood@wustl.edu',
    packages=find_packages(),
    install_requires=[
    ],
    python_requires='>=3.11',
    entry_points={
        'console_scripts': [
            'sortscore=sortscore.run_analysis:main',
        ],
    },
)