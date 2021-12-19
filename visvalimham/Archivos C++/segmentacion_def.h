#include "segmentacion.h"
using namespace std;

Segmentation::Segmentation(string file_name)
: reverse_segment(true), segment_size(0), reduce_size(2), reduce_size_aux(0), jaccard_coef(0)
{   
    build_matrix(file_name);
};

Segmentation::~Segmentation(){
    for (int i=0; i<N1;i++) delete[] image[i];
};

bool Segmentation::is_repited(int y, int x){
    for (vector<Point>::const_iterator rep = segment.begin(); rep != segment.end(); ++rep){
        if ((*rep).first == y && (*rep).second == x) return true;
    }
    return false;
}

void Segmentation::find_extrems(int y, int x){
    if (x < 0 || x >= N2 || y < 0 || y >= N1) return; // no se debe de salir del margen de la matriz
    else if (is_repited(y,x)) return;
    segment.emplace_back(make_pair( y ,x ));
    segment_size++;
    int num_neighbors = 0;
    Table neighbors;
    for (int k = 0; k < 4; k++) { //movimiento en cruz. Importante para evitar bugs en los surcos
        if (image[ y + dy1[k] ][ x+dx1[k] ] > 0){ //esta pintado el pixel?
            if (!is_repited(y+dy1[k], x+dx1[k])){ //no esta repetido en segment?
                num_neighbors++;
                neighbors.emplace_back(make_pair( y+dy1[k], x+dx1[k]));
            }
        }
    }
    if (num_neighbors == 0){
        for (int k =0; k<4; k++){ //movimiento en X si no se encuentra algun vecino en cruz.
            if (image[ y+dy2[k] ][ x+dx2[k] ]>0){
                if (!is_repited(y+dy2[k], x+dx2[k])){
                num_neighbors++;
                neighbors.emplace_back(make_pair(y + dy2[k], x + dx2[k]));
                }
            }
        }
    }
    //si el punto (x,y) no tiene vecinos, es un punto extremo
    if (num_neighbors == 0) fill_segment(); // reacomodamos el vector de segmento. El primer y ultimo elemento son los puntos extremos
    else for(int i =0; i<num_neighbors; i++) find_extrems(neighbors[i].first, neighbors[i].second);
}

void Segmentation::fill_segment(){
    if (reverse_segment){
        reverse(segment.begin(),segment.end());
        reverse_segment = false;
    }
    else return; //en la segunda vuelta, los pixeles se acomodaran en orden en find_extrems
}

void Segmentation::search_segment(){
    //Busca algun punto conexo en toda la imagen
    for (int i = 0; i < N1; i++){
        for (int j=0; j < N2; j++){
            if (image[i][j] > 0) find_extrems(i,j);
        }
    }
}

void Segmentation::build_matrix(string file_name){
  ifstream infile(file_name);
  if (!infile.is_open()) return;
  stringstream ss;
  string inputLine;
  // Primera linea : version P2
  getline(infile,inputLine);
  ss << infile.rdbuf();
  // Segunda linea : dimensiones
  ss >> N1 >> N2;
  // Second line : maximo valor en array
  ss >> max_color;

  image = new int*[N1];
  for (int i =0; i<N1; i++) image[i] = new int[N2];
  
  for(int row = 0; row < N1; ++row){
    for (int col = 0; col < N2; ++col){
        ss >> image[row][col];
    }
  }
  infile.close();
}

double Segmentation::perpendicular_dist(const Point &point, const Point &begin_point, const Point &end_point){
    //Vector paralelo a la recta que une los puntos begin y end
	double dy = (double)begin_point.first - (double)end_point.first;
	double dx = (double)begin_point.second - (double)end_point.second;
    //Normalizando el vector
    double norm = sqrt(pow(dx,2.0)+pow(dy,2.0));
    if(norm > 0.0) dx/=norm, dy/=norm; //condicion para evitar comparar vecinos
    //Vector que apunta al punto point
	double vx = (double)point.second - (double)end_point.second;
	double vy = (double)point.first - (double)end_point.first;
	//Proyeccion del punto pt sobre la recta
	double proy = (dx*vx+dy*vy);
	//Vector perpendicular a la recta y que pasa por point
	vx -= proy * dx;
	vy -= proy * dy;
	return sqrt(pow(vx,2.0)+pow(vy,2.0)); //distancia mas corta a la recta
}

void Segmentation::ramer_douglas_peucker(const Table &inp_segment, double epsilon, Table &out_segment)
{
	if(inp_segment.size()<2) throw invalid_argument("No hay suficientes puntos para la reduccion");
	// Encuentra la diastancia maxima entre el punto begin y end
	double dist_max = 0.0, dist;
	int index_max = 0;
	int end = inp_segment.size()-1;
	for(int i = 1; i < end; i++){
        dist = perpendicular_dist(inp_segment[i], inp_segment[0], inp_segment[end]);
		if (dist > dist_max) index_max = i, dist_max = dist;
	}
	// Si la distancia es mas grande que epsilon, aplica rdp recursivamente
	if(dist_max > epsilon){
		//Divicion de lista
        reduce_size++;
        Table out_slice1;
		Table out_slice2;
		Table inp_slice1(inp_segment.begin(), inp_segment.begin()+index_max+1);
		Table inp_slice2(inp_segment.begin()+index_max, inp_segment.end());
		ramer_douglas_peucker(inp_slice1, epsilon, out_slice1);
		ramer_douglas_peucker(inp_slice2, epsilon, out_slice2);
        
        //Union de listas llenadas
		out_segment.assign(out_slice1.begin(), out_slice1.end()-1);
		out_segment.insert(out_segment.end(), out_slice2.begin(), out_slice2.end());
		if(out_segment.size()<2) throw runtime_error("Surgio algun problema al hacer la reduccion");
	} 
	else 
	{
		//Regresa solo los puntos extremos
		out_segment.clear();
		out_segment.push_back(inp_segment[0]);
		out_segment.push_back(inp_segment[end]);
	}
}

void Segmentation::reduce_with(string method, double epsilon = 5, double delta = 0.5){
    if (method == "rdp") ramer_douglas_peucker(segment, epsilon, reduced_segment);
    else if (method == "visva") visvalingam(epsilon);
    else if (method == "bfls") bfls(segment, reduced_segment, delta);
    else if (method == "bfls_visva") bfls_visvalingam(epsilon, delta);
    else cout<<"Metodo de reduccion '"+ method + "' no esta especificado"<<endl;
}

void Segmentation::save_points(string file_name, char o = 'N'){
    ofstream out_file(file_name);
    if (o == 'x'){
    for(int i=0; i< reduced_segment.size();i++){
		out_file << reduced_segment[i].second << "\t" << -reduced_segment[i].first << endl;
	}}
    else{
        for(int i=0; i< segment.size();i++){
		out_file << segment[i].second << "\t" << -segment[i].first << endl;
    }
    out_file.close();
}

}
void Segmentation::save_image(string file_name, char o = 'N'){
    ofstream out_file(file_name);
    out_file<<"P2"<<endl;
    out_file<<N1<<" "<<N2<<endl;
    out_file<<max_color<<endl;
    if (o =='x'){
        for (int i=0; i<reduced_segment.size(); i++){
            for(int k =0; k < 4; k++){
                for (int sca = 1; sca<3;sca++)
                image[reduced_segment[i].first + sca*dy2[k] ][ reduced_segment[i].second +sca*dx2[k]] = max_color;
            }
        }
    }
    for (int i=0; i<N1; i++) for (int j=0; j<N2; j++) out_file<< image[i][j] <<endl;
    out_file.close();
}

double Segmentation::get_reduc_perc(){
    return 100*(1-(double)reduce_size/segment_size);
}

//-------------------
double Segmentation::parallel_dist(const Point &begin_point, const Point &end_point){
    //Vector paralelo a la recta que une los puntos begin y end
	double dy = (double)begin_point.first - (double)end_point.first;
	double dx = (double)begin_point.second - (double)end_point.second;
    //retorna norma del vector
	return sqrt(pow(dx,2.0)+pow(dy,2.0)); //distancia mas corta entre puntos
}

double Segmentation::triang_area(const Point &point, const Point &begin_point, const Point &end_point){
    return 0.5*perpendicular_dist(point, begin_point, end_point) * parallel_dist(begin_point, end_point); 
}

void Segmentation::visvalingam(double epsilon){
    if(segment_size < 3) throw invalid_argument("No hay suficientes puntos para la reduccion");
	int end = segment_size-1;
    int index_min, index_prev, index_next;
    double min_area;
    double MAX = numeric_limits<double>::max();

    vector<Vertex*> vertex_list(segment_size, NULL); 
    vertex_list[0] = new Vertex(0,0,0,MAX); //primer punto no cuenta
    vertex_list[end] = new Vertex(end,end,end,MAX); //ultimo punto no cuenta

	for(int i = 1; i < end; i++){ //inicializamos los puntos con su area respectiva
        vertex_list[i] =  new Vertex(i,i-1,i+1,triang_area(segment[i], segment[i-1], segment[i+1]));
	}

    do{
	min_area = MAX-1;
    index_min = 0;
    for (int i=1; i<end; i++){
        if(vertex_list[i]->area < min_area){
            min_area = vertex_list[i]->area;
            index_min = i;
        }
    }
    //Eliminamos de la lista index_min
    vertex_list[index_min]->area = MAX;
    index_prev = vertex_list[index_min]->prev;
    index_next = vertex_list[index_min]->next;
    //Actualizamos los vecinos
    vertex_list[index_prev]->next = index_next;
    vertex_list[index_next]->prev = index_prev;
    //Actualizamos el area de los vecinos
    vertex_list[index_prev]->area = triang_area(segment[vertex_list[index_prev]->vertex], segment[vertex_list[index_prev]->prev], segment[vertex_list[index_prev]->next]);      
    vertex_list[index_next]->area = triang_area(segment[vertex_list[index_next]->vertex], segment[vertex_list[index_next]->prev], segment[vertex_list[index_next]->next]);
    } while(min_area < epsilon); //epsilon representa la maxima area deseada
    
    //Agregar los no eliminados
    reduced_segment.emplace_back(segment[0]);
    for(int i=1; i<end; i++){
        if (vertex_list[i]->area != MAX){
        reduced_segment.emplace_back(segment[i]);
        reduce_size++;
        }
    }
    reduced_segment.emplace_back(segment[end]);

    //borrar la memoria heap
    for (int i=0; i<segment_size; i++) delete vertex_list[i];
}
//------------------
double Segmentation::jaccard(){
    int index1 = 0, index2 = 0;
    Point begin, end;
    for (int i=0; i < reduced_segment.size()-1; i++){
        begin = reduced_segment[i]; //puntos a interpolar
        end = reduced_segment[i+1];
        for (int j=index1; j<segment.size()-1; j++){ //buscalos en segement y regresa sus indices
            if (segment[j].first == end.first && segment[j].second==end.second){
                index2 = j;
                break;
            }
        }
        Table slice(segment.begin()+index1, segment.begin()+index2+1);
        jaccard_coef += intersect(slice, begin, end);
        index1=index2;
    }
    return (jaccard_coef-reduce_size+2)/segment_size;
}

double Segmentation::intersect(const Table &list, const Point &begin_point, const Point &end_point){
    double num_inter = 0; //numero de pixeles intersectados
    double dy = (double)begin_point.first - (double)end_point.first;
	double dx = (double)begin_point.second - (double)end_point.second;
    double px, py, y;
    if (dx == 0){
    for (int i=0; i<list.size(); i++){
        py = (double)list[i].first;
        if (dy < 0) if( py-0.5 <= begin_point.first+i && begin_point.first+i <= py+0.5) num_inter++;
        if (dy > 0) if( py-0.5 <= begin_point.first-i && begin_point.first-i <= py+0.5) num_inter++;
    }    
    }
    else{
    for(int i=0; i<list.size(); i++){
        py = (double)list[i].first;
        px = (double)list[i].second;
        y = (dy/dx)*(px-begin_point.second) + begin_point.first;
        if( py-0.5 <= y && y <= py+0.5) num_inter++;
    }
    }
    //cout<<num_inter/list.size()<<" "<<begin_point.second<<" "<<begin_point.first<<endl;
    return num_inter;
}
void Segmentation::bfls_visvalingam(double epsilon, double delta){
    visvalingam(epsilon);
    int index1 = 0, index2=0;
    Point begin, end;
    vector<Table> aux;
    vector<int> indexes;
    for (int i=0; i < reduce_size-1; i++){
        begin = reduced_segment[i]; //puntos a evaluar
        end = reduced_segment[i+1];
        for (int j=index1; j<segment.size(); j++){ //buscalos en segement y regresa sus indices
            if (segment[j].first == end.first && segment[j].second==end.second){
                index2 = j;
                break;
            }
        }
        Table slice(segment.begin()+index1, segment.begin()+index2+1);
        if (intersect(slice, begin, end)/(index2-index1+1) < delta){ //criterio delta > jaccard
            Table new_points;
            bfls(slice, new_points, delta);
            aux.push_back(new_points);
            indexes.push_back(i);
        }
        index1 = index2;
    }
    int desviation = 0;
    for (int i=0; i<indexes.size(); i++){
        reduced_segment.insert(reduced_segment.begin()+indexes[i]+1+desviation, aux[i].begin(), aux[i].end());
        desviation += aux[i].size();
    }
    reduce_size += reduce_size_aux;
}

void Segmentation::bfls(const Table &inp_segment, Table& out_segment, double delta){
    //Best Adjusment Line Segmentation
    double m, b, y, px, py;
    int end = inp_segment.size()-1;
    int mid, index1 = 0, index2 = 0;
    Point point1, point2;
    ///throw invalid_argument("No hay suficientes puntos para la reduccion");
    // Si el coeficiente de jaccard es mayor o igual que delta, regresa los puntos inicial y final.
    // Si no; aplica bals recursivamente
    if(intersect(inp_segment, inp_segment[0], inp_segment[end])/(end+1) >= delta){
        //Regresa solo los puntos extremos
        out_segment.clear();
        out_segment.push_back(inp_segment[0]);
        out_segment.push_back(inp_segment[end]);
        return;
    }
    else if (end < 4){
        mid = end/2;
        out_segment.clear();
        out_segment.push_back(inp_segment[mid]);
        reduce_size_aux--;
        return;
    }
    else{

        regression(inp_segment, m, b); //obtenemos los parametros de bals (m,b)

        index1 = 0;
        for(int i=1; i < end; i++){ //primer punto de interseccion con bals
            py = (double)inp_segment[i].first;
            px = (double)inp_segment[i].second;
            y = m*px + b; //evaluado en la recta de mejor encaje
            if( py-0.5 <= y && y <= py+0.5){
                point1 = inp_segment[i];
                index1 = i;
                break;
            }
        }

        index2 = end;
        for(int i=end-1; i >= 1; i--){ //segundo punto de intersecion con bfls
            py = (double)inp_segment[i].first;
            px = (double)inp_segment[i].second;
            y = m*px + b;
            if( py-0.5 <= y && y <= py+0.5){
                point2 = inp_segment[i];
                index2 = i;
                goto end;
            }
        }

        if (index1 == 0 || index2 == end){ //no se encontro punto de interseccion, retornamos puntos extremos
            out_segment.clear();
            out_segment.push_back(inp_segment[0]);
            out_segment.push_back(inp_segment[end]);
            return;
        }

        end:
        reduce_size_aux +=2;
        Table out_slice1;
        Table out_slice2;
        Table out_slice3;

        //Divicion de lista en 3
        Table inp_slice1(inp_segment.begin(), inp_segment.begin()+index1+1);
        Table inp_slice2(inp_segment.begin()+index1, inp_segment.begin()+index2+1);
        Table inp_slice3(inp_segment.begin()+index2, inp_segment.end());
        //Mandar las listas a bals recursivamente
        bfls(inp_slice1,out_slice1,delta);
        bfls(inp_slice2,out_slice2,delta);
        bfls(inp_slice3,out_slice3,delta);
        
        //Union de listas llenadas
        out_segment.assign(out_slice1.begin(), out_slice1.end()-1);
        out_segment.insert(out_segment.end(), out_slice2.begin(), out_slice2.end()-1);
        out_segment.insert(out_segment.end(), out_slice3.begin(), out_slice3.end());
        if(out_segment.size()<2) throw runtime_error("Surgio algun problema al hacer la reduccion");
    }
}

void Segmentation::regression(const Table &list, double &m, double &b){ //done
    double sumx=0, sumy=0, sumx2=0, sumxy=0, n;
    n = list.size();
    for(int i=0; i < n; i++){
        sumx += (double)list[i].second;
        sumy += (double)list[i].first;
        sumx2 += (double)pow(list[i].second,2.0);
        sumxy += (double)list[i].first * (double)list[i].second;
    }
    m = (n*sumxy - sumx*sumy)/(n*sumx2-pow(sumx,2));
    b = (sumy - m*sumx)/n;
}