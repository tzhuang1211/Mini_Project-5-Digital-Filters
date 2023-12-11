#!/bin/bash

echo "Compiling simple_filter..."
gcc simple_filter.c -o simple_filter -lm || { echo "Compilation failed"; exit 1; }
echo "Compilation successful."

echo "Running simple_filter (8)..."
./simple_filter 8 hL.txt hR.txt YL.txt YR.txt input.wav output.wav || { echo "Execution failed"; exit 1; }
echo "Execution successful."
echo "Running Python script..."
# Use 'python' or 'python3' based on your environment
python show_data.py
echo "Python script execution successful."

echo "Running simple_filter (32)..."
./simple_filter 32 hL.txt hR.txt YL.txt YR.txt input.wav output.wav || { echo "Execution failed"; exit 1; }
echo "Execution successful."
echo "Running Python script..."
# Use 'python' or 'python3' based on your environment
python show_data.py
echo "Python script execution successful."

echo "Running simple_filter (1024)..."
./simple_filter 1024 hL.txt hR.txt YL.txt YR.txt input.wav output.wav || { echo "Execution failed"; exit 1; }
echo "Execution successful."
echo "Running Python script..."
# Use 'python' or 'python3' based on your environment
python show_data.py
echo "Python script execution successful."