

.. _installation:


***************
Installation
***************

Fivepseq is designed to run on Unix systems with python 2.7 or python 3.x. 

Dependencies
------------------

To set up fivepseq, the following python packages need to be pre-installed manually using pip (if you don't have pip you may install it as described `here <https://pip.pypa.io/en/stable/installing/>`_ ). 

Paste the following lines into the shell terminal: 


.. code-block:: shell

 git clone https://github.com/joshuagryphon/plastid -b develop
 cd plastid
 python setup.py install
 pip install --upgrade numpy==1.19.5 pysam==0.19.0 cython==0.29.28


In order to enable exporting vector and portable image files, you'll also need to have phantomjs installed as follows:

.. code-block:: shell

 conda install phantomjs selenium pillow


Fivepseq installation
-----------------------

Clone the project from github: 

.. code-block:: shell

 git clone https://github.com/lilit-nersisyan/fivepseq.git

Navigate into the fivepseq directory and install:

.. code-block:: shell
 
 python setup.py install

To check if fivepseq was installed correctly, type the following in the command line: 

.. code-block:: shell
 
 fivepseq --version

This should display the currently installed version of fivepseq. To display commandline arguments you may type: 

.. code-block:: shell
 
 fivepseq --help


Fivepseq module installation
-----------------------------
In order to install fivepseq as a python module, you may use pip: 

.. code-block:: shell

  pip install fivepseq

