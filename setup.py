import fivepseq
from setuptools import setup, find_packages

def readme():
    with open('README.rst') as f:
        return f.read()

def get_version():
    return fivepseq.__version__

setup(name='fivepseq',
      version=get_version(),
      entry_points={'console_scripts': ['fivepseq = fivepseq.__main__:main']},
      description='A package for analysis and visualization of 5\' endpoint distribution in RNA-seq datasets.',
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
      install_requires=['pathlib2==2.3.7',
                        'preconditions==0.1',
                        'numpy==1.22.0',
                        'pandas==1.1.5',
                        'pysam==0.19.0',
                        'plastid==0.6.1',
                        'dill==0.3.8',
                        'colorlover==0.3.0',
                        'bokeh==2.4.3',
                        'logging==0.4.9.6'],
      zip_safe=False)
