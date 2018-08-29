all:  externalRepos/pbgreedyphase/partitionByPhasedSNVs 

SHELL := /bin/bash

externalRepos/pbgreedyphase/partitionByPhasedSNVs:
	cd externalRepos/pbgreedyphase && make

