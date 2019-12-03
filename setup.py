from setuptools import setup, find_packages

VERSION = '1.0beta0'
import os

def get_version():
    #f = open(os.path.abspath(os.path.join(os.curdir, "fivepseq", "version.txt")), 'r')
    #version = f.read()
    return VERSION

def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='fivepseq',
      version=get_version(),
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
      install_requires=['pathlib2', 'preconditions',
                        'numpy', 'pandas', 'pysam', 'plastid',
                        'dill', 'colorlover',
                        'bokeh', 'logging'],
      zip_safe=False)
