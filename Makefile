all: scripts/readToSNVList \
	envs/python2.done \
	envs/python3.done \
	TestCases/SDAtest/ref.fasta \
	externalRepos/racon-v1.4.5/build/bin/racon \
	externalRepos/canu-1.8/Linux-amd64/bin/canu \

SHELL := /bin/bash

CCOPTS=-g

DEBUG?=""
ifneq ($(DEBUG), "")
	CCOPTS=$(CCOPTS_BASE) $(DEBUG)
else
	CCOPTS=-O3 $(CCOPTS_BASE)
endif
STATIC=-static
ifeq ($(OPT), "1")
	CCOPTS=-g $(CCOPTS_BASE) -lprofiler
	STATIC=
endif



# 
# Get readToSNVList
#
externalRepos/pbgreedyphase/readToSNVList:
	source env_sda.sh && \
	cd externalRepos/pbgreedyphase && \
	make readToSNVList

scripts/readToSNVList: externalRepos/pbgreedyphase/readToSNVList
	source env_sda.sh && \
	cp externalRepos/pbgreedyphase/readToSNVList scripts/readToSNVList

#
# make conda envs
#
envs/python2.done: 
	source env_sda.sh && \
	mkdir -p envs && \
	conda env create --force -f envs/python2.yml --prefix $(PWD)/envs/sda-python-2 && \
	touch envs/python2.done

envs/python3.done:
	source env_sda.sh && \
	mkdir -p envs && \
	conda env create --force -f envs/python3.yml --prefix $(PWD)/envs/sda-python-3 && \
	touch envs/python3.done


#
# make canu 
#
externalRepos/canu-1.8/Linux-amd64/bin/canu:
	source env_sda.sh && \
	cd externalRepos && \
	wget https://github.com/marbl/canu/releases/download/v1.8/canu-1.8.Linux-amd64.tar.xz && \
	tar -xf canu-1.8.Linux-amd64.tar.xz && \
	rm canu-1.8.Linux-amd64.tar.xz

#
# make racon, okay if this fails, the pipeline will run without polishing
#
externalRepos/racon-v1.4.5/build/bin/racon:
	source env_sda.sh && \
	cd externalRepos && \
	wget https://github.com/lbcb-sci/racon/releases/download/1.4.5/racon-v1.4.5.tar.gz && \
	tar -xvzf racon-v1.4.5.tar.gz && \
	cd racon-v1.4.5 && mkdir -p build && cd build && \
	cmake -DCMAKE_BUILD_TYPE=Release .. && make && \
	cd ../../ && rm racon-v1.4.5.tar.gz 


#
# Make hstlib for bamtofreq
#
externalRepos/htslib-1.9/libhts.so:
	source env_sda.sh && \
	cd externalRepos && \
	wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
	tar xjf htslib-1.9.tar.bz2 && \
   	cd htslib-1.9/ && \
	autoheader && autoconf && ./configure --prefix=$(PWD) && make && \
	cd ../ && rm htslib-1.9.tar.bz2

#
#
#
bamToFreq: externalRepos/BamToFreq.cpp externalRepos/htslib-1.9/libhts.so
	source env_sda.sh && \
	g++ -I externalRepos/htslib-1.9 $(CCOPTS) $^ -o $@  -L externalRepos/htslib-1.9 -lhts -lpthread -lz -Wl,-rpath,$(PWD)/externalRepos/htslib-1.9

#
# only make the dependacies required for the collapse analysis. 
#
collapse: envs/python3.done

#
# Make the test cases
#
TestCases/SDAtest/ref.fasta:
	rm -rf TestCases/SDAtest* 
	cd TestCases/ && \
	wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/SDAtest.tar.gz  && \
	wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/SDAtest2.tar.gz  && \
	tar -zxvf SDAtest.tar.gz && \
	tar -zxvf SDAtest2.tar.gz 

#
# Make the fake genome test case
#
TestCases/GenomeTest/ref.fasta:
	rm -rf TestCases/GenomeTest* 
	cd TestCases/ && \
	wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/GenomeTest.tar.gz && \
	tar -zxvf GenomeTest.tar.gz

clean: 
	source env_sda.sh && \
	conda remove -y --prefix $(PWD)/envs/sda-python-3 --all ; \
	conda remove -y --prefix $(PWD)/envs/sda-python-2 --all ; \
	rm -f envs/python2.done envs/python3.done ; \
	rm -rf TestCases/SDAtest* ; \
	rm -rf TestCases/GenomeTest* ; \
	rm -f externalRepos/pbgreedyphase/readToSNVList ; \
	rm -f scripts/readToSNVList


