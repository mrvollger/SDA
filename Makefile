all: setup_abp.sh \
  environments/python2.7/bin/activate \
  quiver/bin/variantCaller.new \
  externalRepos/pbgreedyphase/partitionByPhasedSNVs \
  externalRepos/pbsamstream/pbsamstream \
  environments/python2.7/lib/python2.7/site-packages/h5py-2.8.0.post0-py2.7-linux-x86_64.egg \

SHELL := /bin/bash


externalRepos/htslib/libhts.a:
	cd externalRepos/htslib; \
    autoheader; \
    autoconf; \
    ./configure --disable-bz2 --disable-lzma --disable-libcurl --disable-s3; \
    make -j 4

externalRepos/boost_1_66_0/bootstrap.sh:
	cd externalRepos && wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz && tar xvf boost_1_66_0.tar.gz
	touch $@

externalRepos/boost_1_66_0/stage/lib/libboost_program_options.a: externalRepos/boost_1_66_0/bootstrap.sh
	cd externalRepos/boost_1_66_0 && ./bootstrap.sh --with-libraries=program_options && ./b2 --prefix=$PWD/build -j 4

environments/python2.7/bin/activate:
	./setup_virtualenv.sh


environments/python2.7/lib/python2.7/site-packages/h5py-2.8.0.post0-py2.7-linux-x86_64.egg: environments/python2.7/bin/activate externalRepos/hdf5/build/lib/libhdf5_cpp.so
	source ./environments/python2.7/bin/activate && \
   cd externalRepos/h5py && \
   python setup.py configure --hdf5=$(PWD)/externalRepos/hdf5/build/ && \
   python setup.py build && \
   python setup.py install && \
   python setup.py install_egg_info 


environments/python2.7/lib/python2.7/site-packages/ConsensusCore-1.0.2-py2.7.egg/ConsensusCore.py: externalRepos/boost_1_66_0/stage/lib/libboost_program_options.a environments/python2.7/bin/activate
	mkdir -p quiver/lib/python2.7/site-packages/
	source ./environments/python2.7/bin/activate && \
  cd externalRepos/ConsensusCore && \
  python setup.py build --boost=$(PWD)/externalRepos/boost_1_66_0 && \
  python setup.py install --boost=/home/cmb-16/mjc/mchaisso/projects/abp_v2/externalRepos/boost_1_66_0 


environments/python2.7/lib/python2.7/site-packages/ConsensusCore-1.0.2-py2.7.egg-info:
	source ./environments/python2.7/bin/activate && \
  cd externalRepos/ConsensusCore && python setup.py install_egg_info

environments/python2.7/lib/python2.7/site-packages/pbcommand-1.0.0-py2.7.egg:
	mkdir -p quiver/lib/python2.7/site-packages/
	cd externalRepos/pbcommand && python setup.py build && \
    python setup.py install 

quiver/bin/quiver: environments/python2.7/lib/python2.7/site-packages/ConsensusCore-1.0.2-py2.7.egg/ConsensusCore.py environments/python2.7/lib/python2.7/site-packages/pbcommand-1.0.0-py2.7.egg environments/python2.7/bin/activate environments/python2.7/lib/python2.7/site-packages/h5py-2.8.0.post0-py2.7-linux-x86_64.egg environments/python2.7/lib/python2.7/site-packages/ConsensusCore-1.0.2-py2.7.egg-info
	mkdir -p quiver/lib/python2.7/site-packages/
	source $(PWD)/environments/python2.7/bin/activate && \
    cd externalRepos/GenomicConsensus && \
    python setup.py build && \
    python setup.py install 


externalRepos/hdf5-1.8.14:
	cd externalRepos && wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz && tar zxvf hdf5-1.8.14.tar.gz

externalRepos/hdf5/build/lib/libhdf5_cpp.so: externalRepos/hdf5-1.8.14
	rm -rf externalRepos/hdf5/cmake_build
	mkdir -p externalRepos/hdf5/build
	export CXXFLAGS="-std=c++11"
	cd externalRepos/hdf5-1.8.14 && \
  mkdir cmake_build && \
  cd cmake_build && \
  CXXFLAGS=-std=c++11 && cmake .. -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CPP_COMPILER=`which g++` -DBUILD_SHARED_LIBS:BOOL=ON -DHDF5_BUILD_CPP_LIB:BOOL=ON  -DCMAKE_INSTALL_PREFIX:PATH=$(PWD)/externalRepos/hdf5/build && \
  make -j 8 VERBOSE=1 && \
  make install

externalRepos/pbbam/build/lib/libpbbam.a: externalRepos/hdf5/build/lib/libhdf5_cpp.so externalRepos/boost_1_66_0/bootstrap.sh externalRepos/htslib/libhts.a
	cd externalRepos/pbbam/; \
   rm -rf build; \
   mkdir -p build; cd build; \
   cmake   -D HTSLIB_LIBRARIES=$(PWD)/externalRepos/htslib/libhts.a \
   -DCMAKE_C_COMPILER=`which gcc` \
   -DCMAKE_CPP_COMPILER=`which g++` \
   -D CMAKE_CXX_STANDARD=11 \
   -D HDF5_LIBRARIES=$(PWD)/externalRepos/hdf5/build/lib \
   -D HDF5_INCLUDE_DIRS=$(PWD)/externalRepos/hdf5/build/include \
   -D HTSLIB_INCLUDE_DIRS=$(PWD)/externalRepos/htslib \
   -D PacBioBAM_build_tests=False \
   -D BOOST_ROOT=$(PWD)/externalRepos/boost_1_66_0/ .. ; \
   make VERBOSE=1 -j 8

quiver/bin/pbindex: externalRepos/pbbam/build/lib/libpbbam.a
	mkdir -p quiver/bin
	cp externalRepos/pbbam/build/bin/pbindex quiver/bin/pbindex

quiver/bin/variantCaller.new: quiver/bin/quiver quiver/bin/pbindex
	echo "#!/usr/bin/env python" > quiver/bin/variantCaller.new
	tail -n +2 quiver/bin/variantCaller >> quiver/bin/variantCaller.new
	mv -f quiver/bin/variantCaller.new quiver/bin/variantCaller
	chmod +x  quiver/bin/variantCaller



setup_abp.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "source "$(PWD)"/environments/python2.7/bin/activate" >> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/local_assembly/pbgreedyphase/boost_1_66_0/stage/lib/">> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/quiver/lib/">> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/quiver/lib64/">> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/hdf5/build/lib/">> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/local_assembly/pbgreedyphase/boost_1_66_0/stage/lib/">> $@
	echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:"$(PWD)"/local_assembly/pbgreedyphase/lzma/build/lib">> $@
#	echo "export PYTHONPATH=\$$PYTHONPATH:"$(PWD)"/quiver/lib/python2.7/site-packages/">> $@
	echo "export PATH=\$$PATH:"$(PWD)"/quiver/bin/" >> $@
	echo "export PATH=\$$PATH:"$(PWD)"/bin/" >> $@	
	echo "#" >> $@
	echo "# Add custom configuration here." >> $@
	echo "#" >> $@


hdf5/build/lib/libhdf5_cpp.so:
	rm -rf $(PWD)/hdf5/cmake_build
	mkdir -p $(PWD)/hdf5/build
	export CXXFLAGS="-std=c++11"
	cd hdf5/ && \
  mkdir cmake_build && \
  cd cmake_build && \
  CXXFLAGS=-std=c++11 && cmake .. -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CPP_COMPILER=`which g++` -DBUILD_SHARED_LIBS:BOOL=ON -DHDF5_BUILD_CPP_LIB:BOOL=ON  -DCMAKE_INSTALL_PREFIX:PATH=$(PWD)/hdf5/build && \
  make -j 8 VERBOSE=1 && \
  make install

externalRepos/pbgreedyphase/partitionByPhasedSNVs:
	cd externalRepos/pbgreedyphase && make

externalRepos/pbsamstream/pbsamstream:
	cd externalRepos/pbsamstream && make
