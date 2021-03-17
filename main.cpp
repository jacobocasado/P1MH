// Jacobo Casado de Gracia. Práctica 1 de Metaheurística. Greedy y Búsqueda Local.

/* En primer lugar, la práctica está dividida de estas maneras:
 * 1: Almacenamiento de datos para un mejor uso en los algoritmos.
 * 2: Aplicación del algoritmo Greedy y aplicación del algoritmo BL.
 * 3: Representación y guardado de datos.
 */

using namespace std;

#include <iostream>
#include <fstream>
#include "Eigen/Dense"


Eigen::MatrixXd generarMatrizDistancias(string archivo, int &size){

    int n, m;
    int f, c;
    Eigen::MatrixXd matrizDistancias;

    ifstream lectura;
    lectura.open(archivo, ios::out | ios::in);

    if (lectura.is_open())
    {

        float distancia;

        // Guardamos el tamaño de la matriz y el subconjunto en ambas variables.
        lectura >> n;
        lectura >> m;
        size = m;
        matrizDistancias.resize(n,n);


        while (!lectura.eof()){
            lectura >> f;
            lectura >> c;
            lectura >> distancia;

            matrizDistancias(f,c) = distancia;
            matrizDistancias(c,f) = distancia;
        }
    }
    else
    {
        cout << "El archivo no pudo ser abierto." << endl;
    }
    lectura.close();

    return matrizDistancias;
}

int encontrarPrimerElementoMaximaDistancia(Eigen::MatrixXd matrizDistancias){

    int posicionMejor = -1;
    double distanciaMejor = 0.0;
    double distanciaActual;

    for (int fila = 0; fila < matrizDistancias.rows(); fila++){

        distanciaActual = 0.0;

        for (int col = 0; col < matrizDistancias.cols(); col++){
            distanciaActual += matrizDistancias(fila,col);
        }

        if (distanciaActual > distanciaMejor){
            posicionMejor = fila;
            distanciaMejor = distanciaActual;
        }
    }
    return posicionMejor;
}

int encontrarSiguienteElementoMaximaDistancia(Eigen::MatrixXd matrizDistancias, Eigen::ArrayXi vectorSolucion, int &posicion){

    Eigen::ArrayXi vectorDistancias(matrizDistancias.cols());

    for (int i = 0; i <= posicion; ++i){
        for (int j = 0; j < matrizDistancias.cols(); ++j){
            vectorDistancias(j) += matrizDistancias(i,j);
        }
    }
    // En el vector de distancias estan las distancias acumuladas.
    // Elegimos el maximo del vector, que nos va a decir el siguiente elemento a anadir en cuestion.

    // Obtenemos el maximo
    Eigen::ArrayXi::Index indiceMaximo;
    vectorDistancias.maxCoeff(&indiceMaximo);
    int indice = indiceMaximo;
    return indice;
}


int main() {
    // Lo primero que debemos hacer es obtener los datos de la matriz dada en los archivos de tablas.
    // Probaremos que obtenemos los resultados deseados.

    int size; // Tamanio del subconjunto.

    Eigen::MatrixXd matrizDistancias = generarMatrizDistancias("tablas/MDG-a_1_n500_m50.txt", size);
    Eigen::ArrayXi vectorSolucion(size);

    int primerElemento = encontrarPrimerElementoMaximaDistancia(matrizDistancias);

    vectorSolucion << primerElemento;

    int posicion = 1;

    for (int i = 1; vectorSolucion.size(); i++){
        cout << encontrarSiguienteElementoMaximaDistancia(matrizDistancias, vectorSolucion, posicion);
    }


}
