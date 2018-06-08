all: pbsamstream/hdf5/build/lib/libhdf5.a \
	blasr/alignment/bin/blasr
#all: hdf5/build/lib/libhdf5.so \
	blasr/alignment/bin/blasr



# not used
hdf5/build/lib/libhdf5.so:
	echo "Making HDF5!" && \
		rm -rf $(PWD)/hdf5/build $(PWD)/hdf5/cmake_build && \
		mkdir -p $(PWD)/hdf5/build && \
		export CXXFLAGS="-std=c++11" && \
		cd hdf5/ && \
		mkdir -p cmake_build && \
		cd cmake_build && \
		cmake ..  -DHDF5_BUILD_CPP_LIB:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=ON -DCMAKE_INSTALL_PREFIX:PATH=$(PWD)/hdf5/build  && \
		make -j 8 VERBOSE=1 && \
		make install
		#git checkout -b blasr_hdf5 8a275ab7831002f3e402947facc2d373d9d7926c && 

blasr/alignment/bin/blasr: hdf5/build/lib/libhdf5.so
	echo "Making BLASR!" && \
		export HDF5INCLUDEDIR=$(PWD)/hdf5/build/include/ && \
		export HDF5LIBDIR=$(PWD)/hdf5/build/lib/ && \
		cd blasr && \
		make 


pbsamstream/hdf5/build/lib/libhdf5.a:
	echo "Making pbsamstream/HDF5!" && \
		cd pbsamstream && \
		git submodule update --init --recursive && \
		make && \
		cd ..


blasrx/alignment/bin/blasr: pbsamstream/hdf5/build/lib/libhdf5.a
	echo "Making BLASR!" && \
		export HDF5INCLUDEDIR=$(PWD)/pbsamstream/hdf5/build/include/ && \
		export HDF5LIBDIR=$(PWD)/pbsamstream/hdf5/build/lib/ && \
		cd blasr && \
		make 


