#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>

int Class_Particulas() {
    // Crear un canvas
    TCanvas *c1 = new TCanvas("c1", "Histograma Gaussiano", 800, 600);

    // Crear un histograma de 100 bins entre -5 y 5
    TH1F *h1 = new TH1F("h1", "Distribución Gaussiana;X;Frecuencia", 100, -5, 5);


    // Generador de números aleatorios
    TRandom3 rand(0); // Semilla aleatoria

    // Llenar el histograma con 10000 puntos que siguen una gaussiana de media 0 y sigma 1
    for (int i = 0; i < 10000; ++i) {
        h1->Fill(rand.Gaus(0, 1));
    }

    // Ajustar con una gaussiana
    h1->Fit("gaus");

    // Dibujar el histograma
    h1->Draw();
    
    // Guardar el canvas como PDF
    c1->SaveAs("histograma_gaussiano.pdf");

    return 0;
}
