#!/bin/bash

for f in data/new/*.txt; do
	./bin/wcsp -g $f
done