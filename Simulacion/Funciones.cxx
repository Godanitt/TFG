#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "ActCrossSection.h"
#include "ActTheoCrossSection.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "ActCrossSection.h"
#include "TH1.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"
#include "ROOT/RDataFrame.hxx"

#include <fstream>
#include <vector>
#include <string>
#include <iostream>

void Funciones() {
    // Crear canvas
    TCanvas *c1 = new TCanvas("c1", "Comparacion Breit-Wigner, Gauss y Voigt", 1000, 600);
    
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.05,"XYZ");  // Tamaño de los números en el eje X
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetTitleSize(0.05,"XYZ");  // Tamaño de los números en el eje X
    gStyle->SetTitleFont(62);      // Sin segundo argumento → título general
    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetTitleSize(0.05);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.01);
    gStyle->SetPadTopMargin(0.01);
    gStyle->SetPadBottomMargin(0.12); 
    gROOT->ForceStyle();
    double xmin = -3, xmax = 3;
    double mean = 0.0;
    double sigma = 0.2; // sigma=Gamma

    // Breit-Wigner (Lorentziana)
    TF1 *bw = new TF1("bw", "[0]/(TMath::Pi()*( [1]* (1 + pow((x-[2])/[1], 2))))", xmin, xmax);
    bw->SetParameters(1.0, sigma, mean);
    double integral_bw = bw->Integral(xmin, xmax);
    bw->SetParameter(0, 1.0 / integral_bw);
    bw->SetLineColor(kRed);
    bw->SetLineWidth(3);
    bw->SetNpx(5000); // o más, por defecto es 100


    // Gauss
    TF1 *gauss = new TF1("gauss", "[0]*TMath::Gaus(x, [1], [2], true)", xmin, xmax);
    gauss->SetParameters(1.0, mean, sigma);
    double integral_gauss = gauss->Integral(xmin, xmax);
    gauss->SetParameter(0, 1.0 / integral_gauss);
    gauss->SetLineColor(kBlue);
    gauss->SetLineWidth(3);        
    gauss->SetNpx(5000); // o más, por defecto es 100


    // Voigt (ROOT no permite amplitud en TVoigt, así que definimos manualmente la normalización con una lambda):
    TF1 *voigt = new TF1("voigt", [sigma](double *x, double *p) {
        double val = TMath::Voigt(x[0], sigma, sigma);
        return p[0]*val;
    }, xmin, xmax, 1);
    voigt->SetParameter(0, 1.0); // amplitud inicial
    double integral_voigt = voigt->Integral(xmin, xmax);
    voigt->SetParameter(0, 1.0 / integral_voigt);
    voigt->SetLineColor(kGreen+2);
    voigt->SetLineWidth(3);
    voigt->SetNpx(5000); // o más, por defecto es 100

    // Dibujar marco
    double ymax = std::max(bw->GetMaximum(), std::max(gauss->GetMaximum(), voigt->GetMaximum())) * 1.2;
    TH1D *frame = new TH1D("frame", ";X;f(X)", 100, xmin, xmax);
    frame->SetMinimum(0);
    frame->SetMaximum(ymax);
    frame->Draw();
    

    // Dibujar curvas
    bw->Draw("SAME");
    gauss->Draw("SAME");
    voigt->Draw("SAME");

    // Leyenda
    auto leg = new TLegend(0.7,0.8,0.98,0.98);
    leg->AddEntry(bw, "Breit-Wigner", "l");
    leg->AddEntry(gauss, "Gaussiana", "l");
    leg->AddEntry(voigt, "Voigt", "l");
    leg->Draw();

    gPad->SetLeftMargin(0.2);   // margen izquierdo (fracción del ancho del pad)
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.12);

    // Grosor de los ejes
    gStyle->SetLineWidth(2);      // Más grueso el eje X
    gStyle->SetLineWidth(2);      // Más grueso el eje Y

    // Grosor de los labels
    gStyle->SetLabelSize(0.05,"XYZ");  // Tamaño de los números en el eje X
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetTitleSize(0.05,"XYZ");  // Tamaño de los números en el eje X
    gStyle->SetTitleFont(62);      // Sin segundo argumento → título general
    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetTitleSize(0.05); 

    leg->SetTextSize(0.05);  // Tamaño del texto
    leg->SetTextFont(62);     // Fuente en negrita gStyle->SetPadBottomMargin(0.12);

    // Grosor de los ejes
    gStyle->SetLineWidth(2);      // Más grueso el eje X
    gStyle->SetLineWidth(2);      // Más grueso el eje Y

    // Grosor de los labels

    c1->Update();
    
    c1->SaveAs(TString::Format("/home/daniel/GitHub/TFG/Memoria/Imagenes/PlotVoigtGaussBreitWigner.pdf"));

}
