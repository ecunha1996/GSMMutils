from setuptools import setup, find_packages


with open('README.rst') as readme_file:
    readme = readme_file.read()

if __name__ == "__main__":
    setup(name='gsmmutils',
    version='0.0.3',
    python_requires='>=3.6',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    package_data={"": ["*.xml", "*.csv", "*.txt"]},
    include_package_data=True,
    long_description=readme,
    test_suite='tests'
              )
