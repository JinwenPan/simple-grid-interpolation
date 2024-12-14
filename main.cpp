#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <functional>
#include <sys/time.h>

struct point {
	double x, y, value;//x/y coordinates and the value on this node. 
};

struct triangle {
	int ina, inb, inc;//indices of vertices a/b/c. 
	double bbvax, bbvay, bbvbx, bbvby;//x/y coordinate of bounding box vertex a/b at the bottom-left or top-right corner. 
};

struct StructGrid {
	double originx, originy, meshsize;//the x/y coordinates of the origin and the meshsize h.  
	int cellsnum1d;//the number of cells in 1-dimension, which means the total number of cells is cellsnum1d squared. 
};

struct cell {
	std::vector<int> trind;//a vector of indices of triangles intersecting with this cell. 
};

//functions===============================================================================================================================

//determine whether two rectangles overlap, input x/y coordinates of lower left/upper right vertices of rectangle a and rectangle b. 
//because it's difficult to consider the situations they overlap, just consider the situations they don't overlap. 
//sides overlapping and vertices overlapping are considered as non-overlapping. 
bool RecIns(double alx, double aly, double aux, double auy, double blx, double bly, double bux, double buy) { 
	if ( (blx >= aux) || (bux <= alx) || (bly >= auy) || (buy <= aly) ){	    
		return false;
	}
	else {
		return true;
	}
}

//triangle area formula based on its 3 vertices coordinates
double TriArea(double x1, double y1, double x2, double y2, double x3, double y3) {
	return 0.5 * fabs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
}

//interpolation in a triangle using barycentric coordinates
double InterpolationTri(double px, double py, double x1, double y1, double x2, double y2, double x3, double y3, double v1, double v2, double v3) {
	return ( v1 * TriArea(px, py, x2, y2, x3, y3) / TriArea(x1, y1, x2, y2, x3, y3) +
		     v2 * TriArea(px, py, x3, y3, x1, y1) / TriArea(x1, y1, x2, y2, x3, y3) +
		     v3 * TriArea(px, py, x1, y1, x2, y2) / TriArea(x1, y1, x2, y2, x3, y3) );
}

//Determine whether a point is contained in a triangle.
//vector(p3->p1) cross product with vector(p3->p2), they are in the same plane, so the sign of the result denotes the direction. 
double CrossProduct(double x1, double y1, double x2, double y2, double x3, double y3){
	return (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
}
bool PoInTri(double px, double py, double x1, double y1, double x2, double y2, double x3, double y3){
	bool negative, positive;
	//choose clockwise or anticlockwise direction for a triangle, 
	//if the point is at the same side of all three sides of this triangle, 
	//this point will locate in this triangle. 
	//exist negative cross product.
	negative = (CrossProduct(px, py, x1, y1, x2, y2) < 0) || 
			   (CrossProduct(px, py, x2, y2, x3, y3) < 0) || 
			   (CrossProduct(px, py, x3, y3, x1, y1) < 0)  ;
	//exist positive cross product
	positive = (CrossProduct(px, py, x1, y1, x2, y2) > 0) || 
			   (CrossProduct(px, py, x2, y2, x3, y3) > 0) || 
			   (CrossProduct(px, py, x3, y3, x1, y1) > 0)  ;
	//all these three cross products have the same direction so that the point is in the triangle. 
	return !(negative && positive);
}

//Determine squared distance bwtween two points in 2-dimension. 
double Distance2Points(double x1, double y1, double x2, double y2) {
	return (pow((x1 - x2), 2) + pow((y1 - y2), 2));
}

int main(int argc, char* argv[]) {
	std::string INGRID = argv[1];
	std::string OUTGRID = argv[2];
//read data from ingrid=====================================================
	std::string line;
	std::string word;
//open
	std::ifstream ifsi(INGRID);
	if (!ifsi.is_open()) {
		std::cout << "INGRID File Opening Failure!" << std::endl;
		return 1;
		}
//skip lines
	for (int i = 0; i < 6; i++) {
		getline(ifsi, line);
	}
//read and store the number of points in ingrid. 
	std::istringstream ipointsline(line);
	ipointsline >> word;
	ipointsline >> word;
	int IPointsNum = stoi(word);
//construct a point vector and store x/y coordinates without the data in the third column (2D)
	std::vector<point> IPoints(IPointsNum);
	for (int i = 0; i < IPointsNum; i++) {
		getline(ifsi, line);
		std::istringstream ifsline(line);
		ifsline >> word;
		IPoints[i].x = stod(word);
		ifsline >> word;
		IPoints[i].y = stod(word);
	}
//skip lines and read and store the number or triangles in ingrid. 
	getline(ifsi, line); 
	getline(ifsi, line);
	std::istringstream itrianglesline(line);
	itrianglesline >> word;
	itrianglesline >> word;
	int ITrianglesNum = stoi(word);
//construct a triangle vector and store the indicies of three vertices without storing the unnecessary data of the first column. 
	std::vector<triangle> ITriangles(ITrianglesNum);
	for (int i = 0; i < ITrianglesNum; i++) {
		getline(ifsi, line);
		std::istringstream ifsline(line);
		ifsline >> word;
		ifsline >> word;
		ITriangles[i].ina = stoi(word);
		ifsline >> word;
		ITriangles[i].inb = stoi(word);
		ifsline >> word;
		ITriangles[i].inc = stoi(word);
	}
//skip lines and not store the cell types data. 
	for (int i = 0; i < (ITrianglesNum + 6); i++) {
		getline(ifsi, line);
	}
//read and store the values at each node. 
	for (int i = 0; i < IPointsNum; i++) {
		getline(ifsi, line);
		std::istringstream ifsline(line);
		ifsline >> word;
		IPoints[i].value = stod(word);
	}
//close file
	ifsi.close();

//read data from outgrid========================================================================
//open
	std::ifstream ifso(OUTGRID);
	if (!ifso.is_open()) {
		std::cout << "OURGRID File Opening Failure!" << std::endl;
		return 1;
	}
//skip lines
	for (int i = 0; i < 6; i++) {
		getline(ifso, line);
	}
//read and store the number of points in outgrid. 
	std::istringstream opointsline(line);
	opointsline >> word;
	opointsline >> word;
	int OPointsNum = stoi(word);
//construct a point vector and store x/y coordinates without the data in the third column (2D)
	std::vector<point> OPoints(OPointsNum);
	for (int i = 0; i < OPointsNum; i++) {
		getline(ifso, line);
		std::istringstream ifsline(line);
		ifsline >> word;
		OPoints[i].x = stod(word);
		ifsline >> word;
		OPoints[i].y = stod(word);
	}
//skip lines and read and store the number or triangles in outgrid. 
	getline(ifso, line);
	getline(ifso, line);
	std::istringstream otrianglesline(line);
	otrianglesline >> word;
	otrianglesline >> word;
	int OTrianglesNum = stoi(word);
//construct a triangle vector and store the indicies of three vertices without storing the unnecessary data of the first column. 
	std::vector<triangle> OTriangles(OTrianglesNum);
	for (int i = 0; i < OTrianglesNum; i++) {
		getline(ifso, line);
		std::istringstream ifsline(line);
		ifsline >> word;
		ifsline >> word;
		OTriangles[i].ina = stoi(word);
		ifsline >> word;
		OTriangles[i].inb = stoi(word);
		ifsline >> word;
		OTriangles[i].inc = stoi(word);
	}
//close without storing the values on each node. 
	ifso.close();

//find the largest/smallest x/y coordinates of nodes in ingrid and outgrid========================
	double Minx = IPoints[0].x, Miny = IPoints[0].y, Maxx = IPoints[0].x, Maxy = IPoints[0].y;
//ingrid
	for (int i = 0; i < IPointsNum; i++) {
		if (IPoints[i].x < Minx) {
			Minx = IPoints[i].x;
		}
		if (IPoints[i].y < Miny) {
			Miny = IPoints[i].y;
		}
		if (IPoints[i].x > Maxx) {
			Maxx = IPoints[i].x;
		}
		if (IPoints[i].y > Maxy) {
			Maxy = IPoints[i].y;
		}
	}
//outgrid
	for (int i = 0; i < OPointsNum; i++) {
		if (OPoints[i].x < Minx) {
			Minx = OPoints[i].x;
		}
		if (OPoints[i].y < Miny) {
			Miny = OPoints[i].y;
		}
		if (OPoints[i].x > Maxx) {
			Maxx = OPoints[i].x;
		}
		if (OPoints[i].y > Maxy) {
			Maxy = OPoints[i].y;
		}
	}
//construct a structured buffer grid==========================================================================================================================
	StructGrid BufGrid;
	BufGrid.originx = Minx;//the origin is set at the bottom-left corner. 
	BufGrid.originy = Miny;

	BufGrid.meshsize = 1 / sqrt(IPointsNum);

// start timer --------------------------------
	// int inum = 5;
	// BufGrid.meshsize = pow(2, -inum);
	struct timeval t0, t1;
	gettimeofday(&t0, NULL);
//---------------------------------------------


//choose a longer side to be the side of the structured grid, and plus one is to compensate data type conversion rounding
	BufGrid.cellsnum1d = (((Maxx - Minx) > (Maxy - Miny) ? (Maxx - Minx) : (Maxy - Miny)) / BufGrid.meshsize) + 1;
//construct a bounding box(rectangle) for each triangle in ingrid=============================================================================================
//vertex a is at the bottom-left corner, vertex b is at the top-right corner.
	for (int i = 0; i < ITrianglesNum; i++) {
		ITriangles[i].bbvax = (IPoints[ITriangles[i].ina].x < IPoints[ITriangles[i].inb].x ? IPoints[ITriangles[i].ina].x : IPoints[ITriangles[i].inb].x) 
							  <IPoints[ITriangles[i].inc].x ? 
							  (IPoints[ITriangles[i].ina].x < IPoints[ITriangles[i].inb].x ? IPoints[ITriangles[i].ina].x : IPoints[ITriangles[i].inb].x)
							  :IPoints[ITriangles[i].inc].x;
		ITriangles[i].bbvay = (IPoints[ITriangles[i].ina].y < IPoints[ITriangles[i].inb].y ? IPoints[ITriangles[i].ina].y : IPoints[ITriangles[i].inb].y)
							  <IPoints[ITriangles[i].inc].y ?
							  (IPoints[ITriangles[i].ina].y < IPoints[ITriangles[i].inb].y ? IPoints[ITriangles[i].ina].y : IPoints[ITriangles[i].inb].y)
							  :IPoints[ITriangles[i].inc].y;
		ITriangles[i].bbvbx = (IPoints[ITriangles[i].ina].x > IPoints[ITriangles[i].inb].x ? IPoints[ITriangles[i].ina].x : IPoints[ITriangles[i].inb].x)
							  >IPoints[ITriangles[i].inc].x ?
							  (IPoints[ITriangles[i].ina].x > IPoints[ITriangles[i].inb].x ? IPoints[ITriangles[i].ina].x : IPoints[ITriangles[i].inb].x)
							  :IPoints[ITriangles[i].inc].x;
		ITriangles[i].bbvby = (IPoints[ITriangles[i].ina].y > IPoints[ITriangles[i].inb].y ? IPoints[ITriangles[i].ina].y : IPoints[ITriangles[i].inb].y)
							  >IPoints[ITriangles[i].inc].y ?
							  (IPoints[ITriangles[i].ina].y > IPoints[ITriangles[i].inb].y ? IPoints[ITriangles[i].ina].y : IPoints[ITriangles[i].inb].y)
							  :IPoints[ITriangles[i].inc].y;
	}
//=========================================================================================================
//construct a cell vector to store indices of intersecting triangles for each cell. 
//utilize a pair of ints (i,j)to indicate the location of a cell, 
//which means the most bottom-left cell is denoted as (0,0), 
//the cell at the right of it is (1,0) and the cell at the top of it is (0,1). 

	std::vector<cell> Cell(BufGrid.cellsnum1d * BufGrid.cellsnum1d);
	for (int j = 0; j < BufGrid.cellsnum1d; j++) {
		for (int i = 0; i < BufGrid.cellsnum1d; i++) { //scan each cell
			for (int t = 0; t < ITrianglesNum; t++) { //scan each ingrid triangle and decide whether it overlaps using the bounding box. 
				//if overlapping, push_back the index of the triangle to the vector of this cell. 
				if (RecIns(i * BufGrid.meshsize + BufGrid.originx, j * BufGrid.meshsize + BufGrid.originy, 
					(i + 1) * BufGrid.meshsize + BufGrid.originx, (j + 1) * BufGrid.meshsize + BufGrid.originy,
					ITriangles[t].bbvax, ITriangles[t].bbvay, ITriangles[t].bbvbx, ITriangles[t].bbvby)) {
					Cell[i + j * BufGrid.cellsnum1d].trind.push_back(t);
				}
			}
		}
	}
//interpolation==================================================================================================================
	int i = 0, j = 0;//indices of a cell
	int tri = 0;//index of a triangle
	long unsigned int k = 0;
	for (int m = 0; m < OPointsNum; m++) { //scan each interpolated node and locate it to the cell it belongs to. 
		i = (OPoints[m].x - BufGrid.originx) / BufGrid.meshsize;
		j = (OPoints[m].y - BufGrid.originy) / BufGrid.meshsize;
		//scan each triangle insecting with this cell. 
		for (k = 0; k < Cell[i + j * BufGrid.cellsnum1d].trind.size(); k++) {
			tri = Cell[i + j * BufGrid.cellsnum1d].trind[k];//record the index of a triangle. 
			if (PoInTri(OPoints[m].x, OPoints[m].y,//find the triangle containing the point. 
				IPoints[ITriangles[tri].ina].x, IPoints[ITriangles[tri].ina].y,
				IPoints[ITriangles[tri].inb].x, IPoints[ITriangles[tri].inb].y,
				IPoints[ITriangles[tri].inc].x, IPoints[ITriangles[tri].inc].y)) {
				//interpolation
				OPoints[m].value = InterpolationTri(OPoints[m].x, OPoints[m].y,
				IPoints[ITriangles[tri].ina].x, IPoints[ITriangles[tri].ina].y,
				IPoints[ITriangles[tri].inb].x, IPoints[ITriangles[tri].inb].y,
				IPoints[ITriangles[tri].inc].x, IPoints[ITriangles[tri].inc].y, 
				IPoints[ITriangles[tri].ina].value, 
				IPoints[ITriangles[tri].inb].value,
				IPoints[ITriangles[tri].inc].value );
				break;
			}
		}
		//additionaly, we consider the situation that this point is not contained in all the triangles (bounding box) intersecting with this cell, 
		//which could occur near the boundary of the structured grid. then we need to use another interpolation strategy for this situation: 
		//find the nearest vertices of all the triangles (bounding box) intersecting with this cell to the being interpolated point, 
		//let that vertex's value be the interpolated value of the point. 
		//however, we will ignore the situation that interpolated point is located in a cell which doesn't intersect with any triangle (bounding box). 
		if ((k == Cell[i + j * BufGrid.cellsnum1d].trind.size()) && (Cell[i + j * BufGrid.cellsnum1d].trind.size() != 0)) {
			std::vector<int> PointsIndices; // create a int vector to store the indices of nodes which belong to triangles intersecting with this cell. 
			for (long unsigned int n = 0; n < Cell[i + j * BufGrid.cellsnum1d].trind.size(); n++) {
				tri = Cell[i + j * BufGrid.cellsnum1d].trind[n];
				PointsIndices.push_back(ITriangles[tri].ina);
				PointsIndices.push_back(ITriangles[tri].inb);
				PointsIndices.push_back(ITriangles[tri].inc);
			}
			//delete repeated elements
			sort(PointsIndices.begin(), PointsIndices.end());
			PointsIndices.erase(unique(PointsIndices.begin(), PointsIndices.end()), PointsIndices.end());
			//find the index of nearest points to the interpolated point in this vector. 
			int IndexShortDis = PointsIndices[0];
			for (long unsigned int n = 0; n < PointsIndices.size(); n++) {
				if (Distance2Points(OPoints[m].x, OPoints[m].y, IPoints[PointsIndices[n]].x, IPoints[PointsIndices[n]].y) < 
					Distance2Points(OPoints[m].x, OPoints[m].y, IPoints[IndexShortDis].x, IPoints[IndexShortDis].y)) {
					IndexShortDis = PointsIndices[n];
				}
			}
			OPoints[m].value = IPoints[IndexShortDis].value;
		}
	}

// stop timer --------------------------------
gettimeofday(&t1, NULL);
std::cout << "Wall clock time of interpolation: " <<
((int64_t)(t1.tv_sec - t0.tv_sec) * (int64_t)1000000 +
(int64_t)t1.tv_usec - (int64_t)t0.tv_usec) * 1e-3
<< " ms" << std::endl;
//---------------------------------------------	

//output========================================================================================================================
	std::ofstream ofso("OUTGRID_interpolated.vtk");
	ofso << "# vtk DataFile Version 2.0" << std::endl;
	ofso << "Sample triangulation" << std::endl;
	ofso << "ASCII" << std::endl;
	ofso << "DATASET UNSTRUCTURED_GRID" << std::endl << std::endl;
	ofso << "POINTS " << OPointsNum << " float" << std::endl;
	for (int a = 0; a < OPointsNum; a++) {
		ofso << std::setprecision(12) << OPoints[a].x << " " << OPoints[a].y << " " << "0.0" << std::endl;
	}
	ofso << std::endl << "CELLS " << OTrianglesNum << " " << OTrianglesNum * 4 << std::endl;
	for (int a = 0; a < OTrianglesNum; a++) {
		ofso << 3 << " " << OTriangles[a].ina << " " << OTriangles[a].inb << " " << OTriangles[a].inc << std::endl;
	}
	ofso << std::endl << "CELL_TYPES " << OTrianglesNum << std::endl;
	for (int a = 0; a < OTrianglesNum; a++) {
		ofso << 5 << std::endl;
	}
	ofso << std::endl << "POINT_DATA " << OPointsNum << std::endl;
	ofso << "SCALARS value float 1" << std::endl;
	ofso << "LOOKUP_TABLE default" << std::endl;
	for (int a = 0; a < OPointsNum; a++) {
		ofso << std::setprecision(12) << OPoints[a].value << std::endl;
	}
	ofso << std::endl;
	ofso.close();

	return 0;
}