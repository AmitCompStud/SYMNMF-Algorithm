# SymNMF Clustering Analysis

This project implements Symmetric Non-negative Matrix Factorization (SymNMF) clustering algorithm and compares its performance against K-means clustering using silhouette scores.

## Overview

The project consists of:
- A C extension module (`symnmfmodule`) that provides efficient implementations of SymNMF operations
- Python scripts for running clustering analysis and generating comparative results
- Comparison of SymNMF vs K-means clustering performance on various datasets

## Features

- **SymNMF Operations**: Compute similarity matrix, diagonal degree matrix, normalized similarity matrix, and perform SymNMF factorization
- **Clustering Comparison**: Compare SymNMF and K-means clustering using silhouette scores
- **Batch Analysis**: Generate results for multiple k values and datasets automatically

## Installation

### Prerequisites
- Python 3.x
- GCC compiler (for building the C extension)
- NumPy
- scikit-learn

### Building the C Extension

1. Build the extension module:
```bash
python setup.py build_ext --inplace
```

This will create `symnmfmodule.so` (or `.pyd` on Windows) in the current directory.

### Alternative: Using Makefile

If you want to build the standalone C program:
```bash
make
```

This creates an executable `symnmf` that can be run from command line.

## Usage

### Basic SymNMF Operations

The `symnmf.py` script provides command-line access to SymNMF operations:

```bash
python symnmf.py <k> <goal> <input_file>
```

Where:
- `<k>`: Number of clusters (integer, 2 ≤ k < N where N is number of data points)
- `<goal>`: One of `sym`, `ddg`, `norm`, or `symnmf`
- `<input_file>`: Path to input data file (CSV format, comma-separated floats)

Goals:
- `sym`: Compute similarity matrix A
- `ddg`: Compute diagonal degree matrix D
- `norm`: Compute normalized similarity matrix W
- `symnmf`: Perform full SymNMF factorization and return H matrix

### Clustering Analysis

The `analysis.py` script compares SymNMF and K-means clustering:

```bash
python analysis.py <k> <input_file>
```

Output: Silhouette scores for both algorithms in the format:
```
nmf: 0.XXXX
kmeans: 0.XXXX
```

### Batch Result Generation

The `generate_results.py` script runs analysis for multiple k values and datasets:

```bash
python generate_results.py
```

This will:
- Run analysis for k=2 to 10 on `input_1.txt` and `input_2.txt`
- Run analysis for k=2 to 4 on `input_3.txt`
- Save results to `results.csv`
- Display a summary table

## Input Format

Input files should be CSV files with:
- Each line representing a data point
- Values separated by commas
- All values as floating-point numbers

Example:
```
8.1402,-5.8022,-7.2376
10.1626,-7.4824,-6.5774
9.3153,-5.4974,-6.7025
```

## Output Format

Matrices are output with 4 decimal places, comma-separated values:

```
0.1234,0.5678,0.9012
0.3456,0.7890,0.1234
```

## Project Structure

- `symnmf.c` / `symnmf.h`: Core C implementation of SymNMF algorithms
- `symnmfmodule.c`: Python C extension wrapper
- `symnmf.py`: Python interface and command-line tool
- `kmeans.py`: K-means clustering implementation
- `analysis.py`: Clustering comparison script
- `generate_results.py`: Batch analysis and result generation
- `setup.py`: Build configuration for Python extension
- `Makefile`: Build configuration for standalone C program

## Algorithm Details

### SymNMF Clustering

Symmetric Non-negative Matrix Factorization decomposes a similarity matrix W into H×H^T, where H is a non-negative matrix. The algorithm:

1. Compute similarity matrix A from input data
2. Compute diagonal degree matrix D
3. Compute normalized similarity matrix W = D^(-1/2) × A × D^(-1/2)
4. Initialize H randomly
5. Iteratively update H using multiplicative update rules until convergence

### Silhouette Score

The silhouette score measures how similar an object is to its own cluster vs other clusters. Values range from -1 to 1, where higher values indicate better clustering.

## Dependencies

- numpy
- scikit-learn (for silhouette_score)

Install with:
```bash
pip install numpy scikit-learn
```

## Error Handling

All scripts perform input validation and will exit with "An Error Has Occurred" message for:
- Invalid number of arguments
- Non-numeric k values
- Invalid goal parameters
- k values outside valid range (2 ≤ k < N)
- File reading errors
- Memory allocation failures (in C code)</content>
<parameter name="filePath">vsls:/README.md