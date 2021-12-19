#include <iostream>
#include <string>
#include "segmentacion.h"


using namespace std;

int main(){
    Segmentation imagen_pgm ("lienzo.pgm"); //<-nombre del archivo a segmentar
    imagen_pgm.search_segment();
    //imagen_pgm.save_points("no_reducido.dat");
    imagen_pgm.reduce_with("bfls_visva",3, 0.5); //<-- epsilon delta
    //imagen_pgm.save_points("reducido.dt", 'x');
    //imagen_pgm.save_image("segmentado.pgm", 'x');
    cout<<"El porcentaje de compresion del segmento es "<<imagen_pgm.get_reduc_perc()<<" % "<<endl;
    cout<<"El coeficiente de jaccard es: "<<imagen_pgm.jaccard()<<endl;
	return 0;
}