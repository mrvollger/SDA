export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(echo $(dirname $(which python))/../lib)
virtualenv environments/python2.7 --python=`which python`
source environments/python2.7/bin/activate
p=./environments/python2.7/bin/python
#export PYTHONPATH=$PWD/environments/python2.7/lib/python2.7/site-packages:$PYTHONPATH\
	
$p -m pip install --upgrade pip
$p -m pip install --no-cache-dir six
$p -m pip install --no-cache-dir scipy
$p -m pip install --no-cache-dir h5py
$p -m pip install --no-cache-dir matplotlib
$p -m pip install --no-cache-dir numpy
$p -m pip install --no-cache-dir pysam
$p -m pip install --no-cache-dir networkx
$p -m pip install --no-cache-dir intervaltree
$p -m pip install --no-cache-dir biopython
$p -m pip install --no-cache-dir Cython

