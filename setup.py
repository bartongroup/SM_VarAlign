from distutils.core import setup


setup(
    name='VarAlign',
    version='',
    packages=['varalign', 'varalign.misc'],
    #data_files=[('config', ['varalign/config.txt'])],
    package_data={'varalign': ['config.txt']},
    #include_package_data=True,
    install_requires=['requests', 'pandas', 'pyvcf', 'pysam'],
    url='',
    license='',
    author='smacgowan',
    author_email='s.macgowan@dundee.ac.uk',
    description='This package is used to map and aggregate variants in multiple sequence alignments.'
)
