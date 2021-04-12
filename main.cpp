// Jacobo Casado de Gracia. Práctica 1 de Metaheurística. Greedy y Búsqueda Local.

/* En primer lugar, la práctica está dividida de estas maneras:
 *
 */

using namespace std;

#include <iostream>
#include <fstream>
#include <chrono>
#include <set>
#include <vector>
#include "Eigen/Dense"
#include "random.h"


struct elemento{
    int posicion;
    mutable double diversidad;

    bool operator>(const elemento & elemento) const {

        if (diversidad == elemento.diversidad){
            return (posicion > elemento.posicion);
        }
        else
            return (diversidad > elemento.diversidad);
    }

    // El set va a ordenar con este operador.
    bool operator<(const elemento & elemento) const {

        if (diversidad == elemento.diversidad)
            return posicion < elemento.posicion;
        else
            return (diversidad < elemento.diversidad);
    }


    bool operator==(const elemento & elemento) const {
        return (posicion == elemento.posicion);
    }

    bool operator==(const int & posicion) const {
        return (this->posicion == posicion);
    }

    elemento operator=(const elemento & otro){
        this->diversidad = otro.diversidad;
        this->posicion = otro.posicion;
        return *this;
    }

};

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

double calcularCosteSolucion(vector<elemento> vectorSolucion,Eigen::MatrixXd &matrizDistancias){

    double costeSolucion = 0;

    for (vector<elemento>::iterator it1 = vectorSolucion.begin(); it1 != vectorSolucion.end(); ++it1){
        for (vector<elemento>::iterator it2 = vectorSolucion.begin(); it2 != vectorSolucion.end(); ++it2){
            costeSolucion += matrizDistancias(it1->posicion, it2-> posicion);
        }
    }

    return costeSolucion / 2;
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
    cout << "Tiempo de calculo: " << duration.count() << " segundos" << endl;
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

double calcularContribucionElemento(Eigen::MatrixXd &matrizDistancias, int posicionAIncluir, int posicionAQuitar, vector<elemento> &vectorSolucion){

    double contribucionElemento = 0;

    for (int i = 0; i < vectorSolucion.size(); ++i){
        if (i != posicionAQuitar)
            contribucionElemento += matrizDistancias(posicionAIncluir, vectorSolucion[i].posicion);
    }

    return contribucionElemento;
}

void factorizarSolucion(set<elemento> &setSolucion, Eigen::MatrixXd &matrizDistancias, vector<elemento> &solucionFactorizada){

    solucionFactorizada.resize(setSolucion.size(), {0,0});

    int contador1 = 0;
    int contador2 = 0;

    for (set<elemento>::iterator it1 = setSolucion.begin(); it1 != setSolucion.end(); ++it1){

        contador2 = contador1;
        solucionFactorizada[contador1].posicion = it1->posicion;

        for (set<elemento>::iterator it2 = it1; it2 != setSolucion.end(); ++it2){
            solucionFactorizada[contador1].diversidad += matrizDistancias(it1->posicion, it2->posicion);
            solucionFactorizada[contador2].diversidad += matrizDistancias(it1->posicion, it2->posicion);
            contador2++;
        }

        contador1++;

    }

    sort(solucionFactorizada.begin(), solucionFactorizada.end());

}

void refactorizarVector(Eigen::MatrixXd &matrizDistancias, vector<elemento> &vectorSolucion, int nuevoElemento, int antiguoElemento){

    for (int i = 0; i < vectorSolucion.size(); ++i){
        if(vectorSolucion[i].posicion != nuevoElemento){
            vectorSolucion[i].diversidad -= matrizDistancias(vectorSolucion[i].posicion, antiguoElemento);
            vectorSolucion[i].diversidad += matrizDistancias(vectorSolucion[i].posicion, nuevoElemento);
        }
    }
}

vector<elemento> calcularSolucionBL(Eigen::MatrixXd &matrizDistancias, Eigen::MatrixXd &matrizDistanciasOperadas, int tam){

    set <elemento> setSolucion;
    vector<elemento> vectorSolucion;

    auto start = std::chrono::system_clock::now();

    while (setSolucion.size() < tam){
        elemento elemento;
        elemento.posicion = Randint(0, matrizDistancias.cols() - 1);
        elemento.diversidad = 0;
        setSolucion.insert(elemento);
    }

    // Llenamos el vector de la primera solucion, incluida la factorizacion de todos los elementos.
    factorizarSolucion(setSolucion, matrizDistancias, vectorSolucion );

    int iteraciones = 0;
    const int evaluaciones = 100000;
    bool hayMejoraEnVecindario = true;
    bool hayMejoraIndividual;

    double contribucionElemento;

    vector<elemento>::iterator it = vectorSolucion.begin();

    while (iteraciones < evaluaciones && hayMejoraEnVecindario){

        hayMejoraIndividual = false;

        while (!hayMejoraIndividual && it != vectorSolucion.end())
        {
            for (int i = 0; i < matrizDistancias.cols(); ++i){
                if (find(vectorSolucion.begin(), vectorSolucion.end(), i) == vectorSolucion.end()){
                    contribucionElemento = calcularContribucionElemento(matrizDistancias, i, it->posicion, vectorSolucion);
                    if (contribucionElemento > it->diversidad){
                        elemento elementoMejor = {i, contribucionElemento};
                        hayMejoraIndividual = true;
                        int antiguoElemento = it->posicion;
                        *it = elementoMejor;
                        refactorizarVector(matrizDistancias, vectorSolucion, it->posicion, antiguoElemento);
                    }
                }
            }

            it++;
        }

        if (hayMejoraIndividual){
            sort(vectorSolucion.begin(), vectorSolucion.end());
            it = vectorSolucion.begin();
        }

        if (it == vectorSolucion.end())
            hayMejoraEnVecindario = false;

        iteraciones++;
    }

    auto end = std::chrono::system_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "Tiempo de calculo: " << duration.count() << " segundos" << endl;

    cout << "Coste de la solucion con BL: " << calcularCosteSolucion(vectorSolucion, matrizDistancias);

    return vectorSolucion;

}

int main() {

    // Lo primero que debemos hacer es obtener los datos de la matriz dada en los archivos de tablas.
    // Probaremos que obtenemos los resultados deseados.

    cout.setf(ios::fixed);
    int tam; // Tamanio del subconjunto.

    Eigen::MatrixXd matrizDistancias = generarMatrizDistancias("tablas/MDG-a_2_n500_m50.txt", tam);
    Eigen::MatrixXd matrizDistanciasOperadas = matrizDistancias;
    Eigen::MatrixXd matrizDistanciasOperadas2 = matrizDistancias;

    int semilla = leerDeArchivo("semilla.txt");
    Set_random(semilla);

    // calcularCosteGreedy(matrizDistancias, matrizDistanciasOperadas, tam);
    // calcularCosteGreedy2(matrizDistancias, matrizDistanciasOperadas2, tam);

    calcularSolucionBL(matrizDistancias, matrizDistanciasOperadas, tam);

}
