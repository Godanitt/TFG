/*

En este programa ejecutamos el main, aunque simulation puede ejcutar los archivos de anera independiente, aquí ejecutamos para varias excitaciones.

*/

# include "Simulation.cxx"
#include <fstream>
#include <vector>
#include <string>
#include <iostream>


int main() {

    /*
    Array de índices para seleccionar diferentes tipos de incertidumbre. Cada elemento representa un identificador del tipo de incertidumbre:
    - 0: todas las fuentes de incertidumbre
    - 1: Solo straggling
    - 2: Solo theta
    - 3: ninguna fuente de incertidumbre
    - 4: ninguna fuente de incertidumbre ni briet wigner
    */
    int uncertaintySelector[5] = {0,1,2,3,4};

    std::vector<double> entries = {0.0, 0.0};


    double entriesMatrix[2][2] {0};


    SetMyStyle();
    for (int i = 0; i < 5; i++) {
        std::cout << "Seleccionamos incertidumbre: "  << uncertaintySelector[i] << "\n";
        entries=Simulation(7.5, 0.0, 100000,uncertaintySelector[i]);

        if (i==0){
            entriesMatrix[0][0] = entries[0];
            entriesMatrix[0][1] = entries[1];
            
        }

        entries=Simulation(7.5, 0.20,100000,uncertaintySelector[i]);

        if (i==0){
            entriesMatrix[1][0] = entries[0];
            entriesMatrix[1][1] = entries[1];
            
        }
    }
    
    std::ofstream tex("/home/daniel/GitHub/TFG/Memoria/Cuerpo/TriggerL1.tex");

    tex << "\\begin{tabular}{llll} \\hline\n";  // ajusta columnas según el CSV

    tex << "\\toprule \n";  // ajusta columnas según el CSV
    tex << "$E_{x}=0.0$ MeV & $E_{x}=0.20$ MeV \\\\ \n \\midrule \n";
    std::string linea;
    bool primera = true;
    bool isPrimero=true;
    for (int i=0; i<2; i++) {
        tex<<TString::Format(" %.2f",100.0*entriesMatrix[i][1]/entriesMatrix[i][0])+"\\% &";
    }
    tex << "\\\\ \n\\bottomrule \n";


    tex << "\\end{tabular}\n";



    return 0;
}