#!/bin/bash

# Project: Fish COI DNA Barcoding
# Purpose: Create the project directory structure

echo "Setting up Fish COI project directories..."

mkdir -p data/raw
mkdir -p data/processed
mkdir -p results
mkdir -p logs
mkdir -p docs
mkdir -p scripts

echo "Project directories are ready."

echo "Current project structure:"
find . -maxdepth 2 -type d

