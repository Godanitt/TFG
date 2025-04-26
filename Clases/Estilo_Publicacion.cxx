
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


void SetEstiloPublicacion() {
    gStyle->SetOptStat(0);             // Oculta el recuadro de estadísticas (media, RMS, etc.)
    gStyle->SetOptTitle(0);            // Oculta el título automático del gráfico (puedes usar uno manual)

    gStyle->SetTitleFont(42, "XYZ");   // Fuente moderna (Helvetica-like) para títulos de ejes
    gStyle->SetLabelFont(42, "XYZ");   // Fuente moderna para las etiquetas (números) de los ejes

    gStyle->SetTitleSize(0.05, "XYZ"); // Tamaño del título del eje (más grande que etiquetas)
    gStyle->SetLabelSize(0.045, "XYZ");// Tamaño de las etiquetas del eje (valores numéricos)

    gStyle->SetFrameLineWidth(1);      // Grosor del marco que rodea la gráfica
    gStyle->SetLineWidth(2);           // Grosor por defecto para líneas (funciones, bordes)
    //gStyle->SetHistLineWidth(1.3);       // Grosor de líneas de histogramas (si los usas)
    gStyle->SetTickLength(0.03, "X");  // Aumenta tamaño de los ticks en X
    gStyle->SetTickLength(0.03, "Y");  // Aumenta tamaño de los ticks en Y
    // gStyle->SetEndErrorSize(8);  // tamaño del "sombrerito" del error bar

    gStyle->SetPadLeftMargin(0.15);  // margen izquierdo
    gStyle->SetPadRightMargin(0.03);  // margen derecho
    gStyle->SetPadTopMargin(0.03);    // margen superior
    gStyle->SetPadBottomMargin(0.15); // margen inferior

    gStyle->SetLegendTextSize(0.035); 

}
void SetEstiloPublicacionHisto() {
    gStyle->SetOptStat("nmruoi");             // Oculta el recuadro de estadísticas (media, RMS, etc.)
    //gStyle->SetOptTitle(0);            // Oculta el título automático del gráfico (puedes usar uno manual)

    gStyle->SetTitleFont(42, "XYZ");   // Fuente moderna (Helvetica-like) para títulos de ejes
    gStyle->SetLabelFont(42, "XYZ");   // Fuente moderna para las etiquetas (números) de los ejes

    gStyle->SetTitleSize(0.05, "XYZ"); // Tamaño del título del eje (más grande que etiquetas)
    gStyle->SetLabelSize(0.045, "XYZ");// Tamaño de las etiquetas del eje (valores numéricos)

    gStyle->SetFrameLineWidth(1);      // Grosor del marco que rodea la gráfica
    gStyle->SetLineWidth(2);           // Grosor por defecto para líneas (funciones, bordes)
    //gStyle->SetHistLineWidth(1.3);       // Grosor de líneas de histogramas (si los usas)
    gStyle->SetTickLength(0.03, "X");  // Aumenta tamaño de los ticks en X
    gStyle->SetTickLength(0.03, "Y");  // Aumenta tamaño de los ticks en Y
    // gStyle->SetEndErrorSize(8);  // tamaño del "sombrerito" del error bar

    gStyle->SetPadLeftMargin(0.10);  // margen izquierdo
    gStyle->SetPadRightMargin(0.15);  // margen derecho
    gStyle->SetPadTopMargin(0.03);    // margen superior
    gStyle->SetPadBottomMargin(0.15); // margen inferior

    gStyle->SetLegendTextSize(0.035); 

}

