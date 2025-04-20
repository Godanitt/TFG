#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>


class Particula
{
    double fZ;
    double fA;
    double fmasa_real;

public:
    void definir_AZ(double valorA, double valorZ, double valorM);
    double get_M() { return fmasa_real*931.494; }; // La masa esta multiplicada por el factor de cambio kg --> MeV/c² 
    void Print();
};

void Particula::definir_AZ(double A, double Z, double M)
{
    fA = A;
    fZ = Z;
    fmasa_real = M;
}

void Particula::Print()
{
    std::cout << "--- Particle ----" << '\n';
    std::cout << "-> Z : " << fZ << '\n';
    std::cout << "-> A : " << fA << '\n';
}

