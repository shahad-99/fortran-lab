# Fortran Laboratory Assignments

This repository contains a collection of Fortran programs and subroutines developed as part of a computational methods laboratory course. The code demonstrates various numerical algorithms and data structures implemented in Fortran.

## 📁 Repository Structure

The repository includes the following Fortran programs and modules:

### 🔗 Linked List Programs
- **`Ex-15-1.f95`** - Implementation of a linked list data structure for storing real values
- **`insertion_sort.f90`** - Insertion sort algorithm using linked lists for integer values

### 📊 Array Operations
- **`Ex-15-2.f95`** - Subroutine to extract diagonal elements from 2D arrays with comprehensive error handling
- **`Ex-15-3.f95`** - Customer database management system with multiple sorting options (by last name, city, or zip code)

### 🌳 Binary Tree Implementation
- **`Ex-15-4.f95`** - Binary tree data structure for storing and retrieving personal information (names, phone numbers) with alphabetical sorting

### 🧮 Linear Equation Solvers
- **`Ex9-1.f95`** - Gaussian elimination with maximum pivot technique for solving linear equations
- **`Ex-11-2.f95`** - Enhanced version with double precision arithmetic and automatic arrays to minimize roundoff errors

## 🛠️ Features

### Data Structures
- Dynamic memory allocation with pointers
- Linked lists for sequential data storage
- Binary trees for efficient searching and sorting
- Derived data types for complex data organization

### Numerical Methods
- Gaussian elimination with partial pivoting
- Insertion sort algorithm
- Matrix operations and diagonal extraction
- Single and double precision numerical solvers

### Error Handling
- Comprehensive error checking
- Memory allocation status verification
- Input validation and file operation checks

## 📋 Requirements

- Fortran 90/95 compiler (gfortran, ifort, etc.)
- Basic understanding of Fortran programming
- Input data files for programs that require external data

## 🚀 Compilation and Execution

### Basic Compilation
```bash
gfortran -o program_name program_name.f90
./program_name
