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


void ponerACeroFila(Eigen::MatrixXd& matriz, unsigned int numFilaARemover)
{
    unsigned int numFilas = matriz.rows();

    if( numFilaARemover < numFilas )
        for (int fila = 0; fila < matriz.rows(); ++fila)
            matriz(numFilaARemover, fila) = 0;

}


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

int encontrarPrimerElementoMaximaDistancia(Eigen::MatrixXd &matrizDistancias){

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
    ponerACeroFila(matrizDistancias, posicionMejor);

    return posicionMejor;
}

void encontrarSiguienteElementoMaximaDistancia(Eigen::MatrixXd &matrizDistancias, Eigen::ArrayXi &vectorSolucion, int &aniadidos){

    Eigen::ArrayXi vectorDistancias(matrizDistancias.rows());
    vectorDistancias.fill(0);


    for (int i = 0; i < aniadidos; ++i){
        for (int j = 0; j < matrizDistancias.rows(); ++j){
            vectorDistancias(j) += matrizDistancias(j,vectorSolucion(i));
        }
    }

    // En el vector de distancias estan las distancias acumuladas.
    // Elegimos el maximo del vector, que nos va a decir el siguiente elemento a anadir en cuestion.

    // Obtenemos el maximo
    Eigen::ArrayXi::Index indiceMaximo;
    vectorDistancias.maxCoeff(&indiceMaximo);
    int indice = indiceMaximo;
    ponerACeroFila(matrizDistancias, indice);

    vectorSolucion(aniadidos) = indice;
    aniadidos++;


}

double calcularCosteTotal(Eigen::ArrayXi vectorSolucion,Eigen::MatrixXd matrizDistancias){

    double distanciaTotal = 0;

    for (int i = 0; i < vectorSolucion.size(); ++i){
        for (int j = i+1; j < vectorSolucion.size(); ++j){
            distanciaTotal += matrizDistancias(vectorSolucion(i),vectorSolucion(j));
        }
    }

    return distanciaTotal;
}


int main() {
    // Lo primero que debemos hacer es obtener los datos de la matriz dada en los archivos de tablas.
    // Probaremos que obtenemos los resultados deseados.

    int tam; // Tamanio del subconjunto.

    Eigen::MatrixXd matrizDistancias = generarMatrizDistancias("tablas/MDG-b_30_n2000_m200.txt", tam);
    Eigen::MatrixXd matrizDistanciasOperadas = matrizDistancias;
    Eigen::ArrayXi vectorSolucion(tam);


    int primerElemento = encontrarPrimerElementoMaximaDistancia(matrizDistanciasOperadas);

    vectorSolucion.fill(0);
    vectorSolucion(0)  = primerElemento;

    int aniadidos = 1;

    while (aniadidos < tam)
        encontrarSiguienteElementoMaximaDistancia(matrizDistanciasOperadas, vectorSolucion, aniadidos);

    cout << calcularCosteTotal(vectorSolucion, matrizDistancias);


}
