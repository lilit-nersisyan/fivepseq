from setuptools import setup, find_packages

VERSION = '0.1.6'


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='fivepseq',
      version=VERSION,
      entry_points={'console_scripts': ['fivepseq = fivepseq.__main__:main']},
      description='A package for analysis of 5pseq datasets',
      url='http://github.com/lilit-nersisyan/fivepseq',
      dependency_links=['http://github.com/lilit-nersisyan/fivepseq/tarball/master#egg=package-1.0'],
      author='Lilit Nersisyan, Maryia Ropat',
      author_email='lilit.nersisyan@ki.se',
      classifiers=[
          'Development Status :: 1 - Planning',
          'Environment :: Win32 (MS Windows)',
          'Intended Audience :: Science/Research',
          'License :: Free for non-commercial use'
          'License :: OSI Approved :: BSD License',
          'Natural Language :: English',
          'Operating System :: Unix',
          'Operating System :: Microsoft :: Windows :: Windows 10',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      license='BSD 3-Clause',
      packages=find_packages(),
      include_package_data=True,
      test_suite='nose.collector',
      tests_require=['nose'],
      install_requires=['pathlib2', 'preconditions', 'plastid', 'numpy', 'pandas', 'pysam', 'dill', 'colorlover',
                        'bokeh', 'logging'],
      zip_safe=False)
