// Jacobo Casado de Gracia. Práctica 1 de Metaheurística. Greedy y Búsqueda Local.

/* En primer lugar, la práctica está dividida de estas maneras:
 *
 */

using namespace std;

#include <iostream>
#include <fstream>
#include <chrono>
#include "Eigen/Dense"
#include "random.h"
#include <set>


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

int leerDeArchivo(string archivo){

    ifstream lectura;
    lectura.open(archivo, ios::out | ios::in);
    int semilla;

    if (lectura.is_open())
    {
        // Guardamos la semilla en la variable.
        lectura >> semilla;
    }
    else
    {
        cout << "El archivo no pudo ser abierto." << endl;
    }
    lectura.close();

    return semilla;
}

void ponerACeroFila(Eigen::MatrixXd& matriz, unsigned int numFilaARemover)
{
    unsigned int numFilas = matriz.rows();

    if( numFilaARemover < numFilas )
        for (int fila = 0; fila < matriz.rows(); ++fila)
            matriz(numFilaARemover, fila) = 0;

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

    int posicionMejor = -1;
    double distanciaMejor = 0.0;
    double distanciaActual;



        for (int j = 0; j < matrizDistancias.rows(); ++j){
            distanciaActual = 0;
            for (int i = 0; i < aniadidos; ++i){
                distanciaActual += matrizDistancias(j, vectorSolucion(i));
            }
             if (distanciaActual > distanciaMejor){
                 posicionMejor = j;
                 distanciaMejor = distanciaActual;
             }
        }

    ponerACeroFila(matrizDistancias, posicionMejor);

    vectorSolucion(aniadidos) = posicionMejor;
    aniadidos++;

}

void encontrarSiguienteElementoMaximaDistancia2(Eigen::MatrixXd &matrizDistancias, Eigen::ArrayXi &vectorSolucion, int &aniadidos){

    Eigen::ArrayXi vectorDistancias(matrizDistancias.rows());
    vectorDistancias.fill(0);

    int posicionMejor = -1;
    double distanciaMejor = 0.0;
    double distanciaMinimaElementoActual;


    for (int i = 0; i < matrizDistancias.rows(); ++i){
        distanciaMinimaElementoActual = matrizDistancias(i, vectorSolucion(0));
        for (int j = 0; j < aniadidos; ++j){
            if (matrizDistancias(i, vectorSolucion(j)) < distanciaMinimaElementoActual){
                distanciaMinimaElementoActual = matrizDistancias(i, vectorSolucion(j));
            }
        }

        if (distanciaMinimaElementoActual > distanciaMejor){
            posicionMejor = i;
            distanciaMejor = distanciaMinimaElementoActual;
        }
    }

    ponerACeroFila(matrizDistancias, posicionMejor);

    vectorSolucion(aniadidos) = posicionMejor;
    aniadidos++;

}

double calcularCosteTotal(Eigen::ArrayXi vectorSolucion,Eigen::MatrixXd &matrizDistancias){

    double distanciaTotal = 0.0;

    for (int i = 0; i < vectorSolucion.size(); ++i){
        for (int j = i+1; j < vectorSolucion.size(); ++j){
            distanciaTotal += matrizDistancias(vectorSolucion(i),vectorSolucion(j));
        }
    }

    return distanciaTotal;
}

double calcularCosteGreedy(Eigen::MatrixXd &matrizDistancias, Eigen::MatrixXd &matrizDistanciasOperadas, int tam ){

    Eigen::ArrayXi vectorSolucion(tam);
    vectorSolucion.fill(0);

    auto start = std::chrono::system_clock::now();

    int primerElemento = encontrarPrimerElementoMaximaDistancia(matrizDistanciasOperadas);
    vectorSolucion(0)  = primerElemento;

    int aniadidos = 1;

    while (aniadidos < tam)
        encontrarSiguienteElementoMaximaDistancia(matrizDistanciasOperadas, vectorSolucion, aniadidos);

    double costeTotalGreedy = calcularCosteTotal(vectorSolucion, matrizDistancias);

    auto end = std::chrono::system_clock::now();
    chrono::duration<double> duration = end - start;

    cout << "Coste Total con Greedy: " << costeTotalGreedy << endl;
    cout << "Tiempo de cálculo: " << duration.count() << " segundos" << endl;
}

double calcularCosteGreedy2(Eigen::MatrixXd &matrizDistancias, Eigen::MatrixXd &matrizDistanciasOperadas, int tam ){

    Eigen::ArrayXi vectorSolucion(tam);
    vectorSolucion.fill(0);

    auto start = std::chrono::system_clock::now();

    int primerElemento = encontrarPrimerElementoMaximaDistancia(matrizDistanciasOperadas);
    vectorSolucion(0)  = primerElemento;

    int aniadidos = 1;

    while (aniadidos < tam)
        encontrarSiguienteElementoMaximaDistancia2(matrizDistanciasOperadas, vectorSolucion, aniadidos);

    double costeTotalGreedy = calcularCosteTotal(vectorSolucion, matrizDistancias);

    auto end = std::chrono::system_clock::now();
    chrono::duration<double> duration = end - start;

    cout << "Coste Total con Greedy: " << costeTotalGreedy << endl;
    cout << "Tiempo de cálculo: " << duration.count() << " segundos" << endl;
}

double calcularCosteBusquedaLocal(Eigen::MatrixXd &matrizDistancias, Eigen::MatrixXd &matrizDistanciasOperadas, int tam){

    struct elemento{
        int posicion;
        double diversidad;

        bool operator>(const elemento & elemento) const {

            if (diversidad == elemento.diversidad){
                if (posicion > elemento.posicion)
                    return false;
                else
                    return true;
            }
            else
                return (diversidad > elemento.diversidad);
        }

        // El set va a ordenar con este operador.
        bool operator<(const elemento & elemento) const {
            if (posicion == elemento.posicion)
                return false;

            if (diversidad == elemento.diversidad)
                return posicion < elemento.posicion;
            else
                return (diversidad < elemento.diversidad);
        }

        bool operator=(const elemento & elemento) const {
            return (posicion == elemento.posicion);
        }

    };

    set <elemento> setSolucion;

    while (setSolucion.size() < tam){
        elemento elemento;
        elemento.posicion = Randint(0, matrizDistancias.cols() - 1);
        elemento.diversidad = 0;
        setSolucion.insert(elemento);
        cout << "Inserto. " << setSolucion.size() << endl;
    }

    for (set<elemento>::iterator it = setSolucion.begin(); it != setSolucion.end(); ++it){
        cout << (*it).posicion << endl;
    }
}

double factorizarSolucion();

int main() {
    // Lo primero que debemos hacer es obtener los datos de la matriz dada en los archivos de tablas.
    // Probaremos que obtenemos los resultados deseados.

    cout.setf(ios::fixed);
    int tam; // Tamanio del subconjunto.

    Eigen::MatrixXd matrizDistancias = generarMatrizDistancias("tablas/MDG-a_1_n500_m50.txt", tam);
    Eigen::MatrixXd matrizDistanciasOperadas = matrizDistancias;
    Eigen::MatrixXd matrizDistanciasOperadas2 = matrizDistancias;

    int semilla = leerDeArchivo("semilla.txt");
    Set_random(semilla);

    // calcularCosteGreedy(matrizDistancias, matrizDistanciasOperadas, tam);
    // calcularCosteGreedy2(matrizDistancias, matrizDistanciasOperadas2, tam);

    calcularCosteBusquedaLocal(matrizDistancias, matrizDistanciasOperadas, tam);

}
