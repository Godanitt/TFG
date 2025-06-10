/*

En este programa ejecutamos el main, aunque simulation puede ejcutar los archivos de anera independiente, aquí ejecutamos para varias excitaciones.

*/

# include "Simulation.cxx"


int main() {

    /*
    Array de índices para seleccionar diferentes tipos de incertidumbre. Cada elemento representa un identificador del tipo de incertidumbre:
    - 0: todas las fuentes de incertidumbre
    - 1: Solo straggling
    - 2: Solo theta
    - 3: ninguna fuente de incertidumbre
    */
    int uncertaintySelector[4] = {0,1,2,3};

    for (int i = 0; i < 4; i++) {
        std::cout << "Seleccionamos incertidumbre: "  << uncertaintySelector[i] << "\n";
        Simulation(7.5, 0.0, 100000,uncertaintySelector[i]);
        Simulation(7.5, 0.20,100000,uncertaintySelector[i]);
    }



    return 0;
}