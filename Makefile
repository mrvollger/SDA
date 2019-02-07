all: scripts/readToSNVList \
	envs/python2.done \
	envs/python3.done \
	TestCases/SDAtest/ref.fasta

SHELL := /bin/bash


# 
# Get readToSNVList
#
externalRepos/pbgreedyphase/readToSNVList:
	source env_conda.cfg && \
		cd externalRepos/pbgreedyphase && \
		make readToSNVList

scripts/readToSNVList: externalRepos/pbgreedyphase/readToSNVList
	source env_conda.cfg && \
		cp externalRepos/pbgreedyphase/readToSNVList scripts/readToSNVList

#
# make conda envs
#
envs/python2.done: 
	source env_conda.cfg && \
		mkdir -p envs && \
	   	conda env create --force -f envs/python2.yml --prefix $(PWD)/envs/sda-python-2 && \
	   	touch envs/python2.done

envs/python3.done:
	source env_conda.cfg && \
		mkdir -p envs && \
	   	conda env create --force -f envs/python3.yml --prefix $(PWD)/envs/sda-python-3 && \
	   	touch envs/python3.done

#
# Make the test case
#
TestCases/SDAtest/ref.fasta:
	rm -rf TestCases/SDAtest* 
	cd TestCases/ && \
	wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/SDAtest.tar.gz  && \
	wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/SDAtest2.tar.gz  && \
	tar -zxvf SDAtest.tar.gz && \
	tar -zxvf SDAtest2.tar.gz 

clean: 
	source env_conda.cfg && \
	conda remove -y --prefix $(PWD)/envs/sda-python-3 --all ; \
	conda remove -y --prefix $(PWD)/envs/sda-python-2 --all ; \
	rm -f envs/python2.done envs/python3.done ; \
	rm -rf TestCases/SDAtest* ; \
	rm -f externalRepos/pbgreedyphase/readToSNVList ; \
	rm -f scripts/readToSNVList

