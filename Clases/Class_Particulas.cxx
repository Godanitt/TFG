
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPaveStats.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Clase que define a una partícula, determinado por su número atómico (Z), número de nucleones (A), masa real (M) y energía de excitación (E) (de momento).
class Particula
{
    double fZ;
    double fA;
    double fmasa_real;
    double fexc;

public:
    // Constructor con parámetros
    void definir_AZ(double valorA, double valorZ, double valorM, double valorExc);

     // La masa está en MeV/c² 
    double get_M() { return fmasa_real; };

     // Z es el número atómico
    double get_Z() { return fZ; };

     // A es el número de nucleones
    double get_A() { return fA; };
     // Exc es la energía de excitación
    double get_Exc() { return fexc; };


    void Print();
    
    // Para poder igualar una partícula a otra
    Particula& operator=(const Particula& other) = default;
};

// Definimos la función definir_AZ() para asignar valores a las variables
// Esta función recibe cuatro parámetros:
// - A (número de nucleones),
// - Z (número atómico), 
// - DM (exceso de masa) [MeV] 
// - E (energía de excitación) [MeV]
// La función asigna estos valores a las variables miembro de la clase Particula
void Particula::definir_AZ(double A, double Z, double DM, double E)
{
    fA = A;
    fZ = Z;
    fmasa_real = A*931.494+DM+E; // La masa está multiplicada por el factor de cambio  u --> MeV/c²
    fexc = E;
}

// Definimos la función Print() para mostrar los valores de las variables
void Particula::Print()
{
    std::cout << "--- Particle ----" << '\n';
    std::cout << "-> Z : " << fZ << '\n';
    std::cout << "-> A : " << fA << '\n';
    std::cout << "-> Ex : " << fexc << '\n';
}

