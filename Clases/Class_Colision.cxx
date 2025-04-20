
#include "Class_Particulas.cxx"

///////////////////////////////////////////////
/*  

Una colisión esta caracterizada por la energía del atomo incidente y las partículas que la componen, por lo que es necesario que formen parte de la clase colisión. Otras que caracterizarían la clase colisión serían: la energía incidente, los ángulos de salida/energías de las partículas resultantes. 

*/
///////////////////////////////////////////////

class Colision
{
    Particula fp1; 
    Particula fp2;
    Particula fp3;
    Particula fp4;
    double ftbeam;

public:
    void set_energia(double valor);
    void set_particulas(Particula p1, Particula p2, Particula p3, Particula p4);

    std::pair<double, double> valores3(Particula p1, Particula p2, Particula p3, Particula p4, double tbeam1, double thetap);
    std::pair<double, double> valores4(Particula p1, Particula p2, Particula p3, Particula p4, double tbeam1, double thetap);

    double get_m1() { return fp1.get_M(); };
    double get_m2() { return fp2.get_M(); };
    double get_m3() { return fp3.get_M(); };
    double get_m4() { return fp4.get_M(); };
    double get_tbeam() { return ftbeam; };
    Particula GetParticle1() { return fp1; }
};


void Colision::set_particulas(Particula p1, Particula p2, Particula p3, Particula p4)
{
    fp1 = p1;
    fp2 = p2;
    fp3 = p3;
    fp4 = p4;
};

void Colision::set_energia(double valor)
{
    ftbeam = valor;
};


std::pair<double, double> valores4(Particula par1, Particula par2, Particula par3, Particula par4, double tbeam1, double thetap){

    double m1 {par1.get_M()};
    double m2 {par2.get_M()};
    double m3 {par3.get_M()};
    double m4 {par4.get_M()};
    double t1 {tbeam1};

    double E1 {m1 + t1};
    double p1 {pow(pow(E1, 2) - pow(m1, 2), 0.5)};
    double E2 {m2};
    double Etot {E1 + m2};
    double Etotp {pow(pow(Etot, 2) - pow(p1, 2), 0.5)};

    double gamma {Etot / Etotp};
    double beta {p1 / Etot};


    double E3p {(Etotp + (pow(m3, 2) - pow(m4, 2)) / Etotp) / 2};
    double E4p {(Etotp + (pow(m4, 2) - pow(m3, 2)) / Etotp) / 2};

    double p4p {pow(pow(E4p, 2) - pow(m4, 2), 0.5)};

    double E4 {gamma * (E4p + beta * p4p * TMath::Cos(thetap))};
    double p4 {pow(-pow(m4, 2) + pow(E4, 2), 0.5)};

    double T4 {E4 - m4};
    double theta4 {TMath::ACos((gamma / p4) * (p4p * TMath::Cos(thetap) + beta * E4p))};
    return {T4,theta4};
}

std::pair<double, double> valores3(Particula par1, Particula par2, Particula par3, Particula par4, double tbeam1, double thetap){

    double m1 {par1.get_M()};
    double m2 {par2.get_M()};
    double m3 {par3.get_M()};
    double m4 {par4.get_M()};
    double t1 {tbeam1};

    double E1 {m1 + t1};
    double p1 {pow(pow(E1, 2) - pow(m1, 2), 0.5)};
    double E2 {m2};
    double Etot {E1 + m2};
    double Etotp {pow(pow(Etot, 2) - pow(p1, 2), 0.5)};

    double gamma {Etot / Etotp};
    double beta {p1 / Etot};


    double E3p {(Etotp + (pow(m3, 2) - pow(m4, 2)) / Etotp) / 2};
    double E4p {(Etotp + (pow(m4, 2) - pow(m3, 2)) / Etotp) / 2};

    double p3p {pow(pow(E3p, 2) - pow(m3, 2), 0.5)};

    double E3 {gamma * (E3p + beta * p3p * TMath::Cos(thetap))};
    double p3 {pow(-pow(m3, 2) + pow(E3, 2), 0.5)};

    double T3 {E3 - m3};
    double theta3 {TMath::ACos((gamma / p3) * (p3p * TMath::Cos(thetap) + beta * E3p))};
    return {T3,theta3};
}

