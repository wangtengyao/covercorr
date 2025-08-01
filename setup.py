from setuptools import setup, find_packages

setup(
    name='covercorr',
    version='1.0.0',
    description='Coverage correlation coefficient implementation for two-sample independence testing',
    author='Tengyao Wang',
    author_email='t.wang59@lse.ac.uk',
    url='https://github.com/wangtengyao/covercorr',
    package_dir={"": "python"},
    packages=find_packages(where="python"),
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'matplotlib',
    ],
    python_requires='>=3.7',
)

