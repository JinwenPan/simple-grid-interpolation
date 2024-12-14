# Grid Interpolation

This repository provides a code implementation for **grid interpolation** between unstructured grids. The code supports interpolation from one unstructured grid with data to another unstructured grid without data. The code processes input and output grids in VTK format and generates interpolated results.

## Requirements
- A C++ compiler supporting the `make` build system.
- [ParaView](https://www.paraview.org/) or similar software to visualize VTK files.

## Getting Started

### 1. Setup
1. Clone or download this repository.
2. Place all files in the same directory.

### 2. Build the Program
Run the following command in the terminal:
```bash
make
```
This will generate an executable file named `intug` in the same directory.

### 3. Run the Program
To perform grid interpolation, execute:
```bash
./intug mesh*in.vtk mesh*out.vtk
```

where:
- `mesh*in.vtk` is the input mesh file containing data.
- `mesh*out.vtk` is the target mesh file without data.

The output will be saved in a new file named `OUTGRID_interpolated.vtk`. All the `vtk` files can be visualized with Paraview. 

### 4. Example Files
- The repository includes several example mesh files named `mesh*.vtk`.
- These files contain the mesh data (for input) and the target mesh without data (for interpolation).

### 5. Additional Resources
- **Slides**: There is a set of slides included in this repository explaining the concept of grid interpolation. The slides provide detailed information on how the interpolation is performed and the underlying methodology.

## Notes
- Ensure all required files are in the same directory for proper execution.
- The program is designed for simplicity and educational purposes. Feel free to modify or extend it for your specific needs.

---

Enjoy exploring grid interpolation with this lightweight tool!

