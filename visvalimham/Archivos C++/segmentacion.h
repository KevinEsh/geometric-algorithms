#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <limits>
using namespace std;

typedef pair<int, int> Point;
typedef vector<Point> Table;

const int dy1[] = {1,-1,0,0};  //movimiento en cruz
const int dx1[] = {0,0,1,-1};  //vecinos del par (x,y): N,S,E,O .Importante este orden
const int dy2[] = {1,-1,-1,1}; //movimiento en X.
const int dx2[] = {1,1,-1,-1}; //vecinos del par (x,y): NE,SE,SO,NO

struct Vertex
{
    Vertex(int vertex_, int prev_, int next_, double area_)
        : vertex(vertex_), prev(prev_), next(next_), area(area_)
    {
    }
    int vertex, prev, next; //un triangulo
    double area; // area del triangulo
};

class Segmentation
{
private:
    int N1, N2, max_color; //dimension de la imagen FILAS, COLUMNAS
    bool reverse_segment;
    int segment_size, reduce_size, reduce_size_aux;
    double jaccard_coef;
    int** image;
    Table segment, reduced_segment;
    void build_matrix(string);
    void find_extrems(int, int);
    void fill_segment();
    bool is_repited(int, int);
    double triang_area(const Point&, const Point&, const Point&);
    double perpendicular_dist(const Point&, const Point&, const Point&);
    double parallel_dist(const Point&, const Point&);
    void ramer_douglas_peucker(const Table&, double, Table&);
    void visvalingam(double);
    void bfls(const Table&, Table&, double);
    void bfls_visvalingam(double, double);
    void regression(const Table&, double&, double&);
    double intersect(const Table&, const Point&, const Point&);
public:
    Segmentation(string);
    ~Segmentation();
    void reduce_with(string, double, double);
    void save_points(string, char);
    void search_segment(); //modificarlo para el proyecto parcial
    void save_image(string, char);
    double reduction_percentage();
    double get_reduc_perc();
    double jaccard();
};

#include "segmentacion_def.h"

#endif