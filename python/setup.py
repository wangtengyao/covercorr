from setuptools import setup, find_packages

setup(
    name='covercorr',
    version='1.0.0',
    description='Coverage correlation coefficient implementation for two-sample independence testing',
    author='Tengyao Wang',
    author_email='t.wang59@lse.ac.uk',
    url='https://github.com/wangtengyao/covercorr',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'scipy',
    ],
    python_requires='>=3.7',
)

