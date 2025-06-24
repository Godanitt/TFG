#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

void ExportSigmasToCSV(const double sigma[2][4],const double usigma[2][4], const std::string& filename = "sigmas.csv") {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "No se pudo abrir el archivo para escritura.\n";
        return;
    }

    // Cabecera (opcional)
    file << ",$\\sigma(0.0)$ [MeV] , $\\sigma(0.20)$ [MeV] \n";

    std::string* firstColumn[4];  // array de 4 punteros a std::string

    firstColumn[0] = new std::string("$\\sigma_{tot}$");
    firstColumn[1] = new std::string("$\\sigma_{straggling}$");
    firstColumn[2] = new std::string("$\\sigma_{\\theta}$");
    firstColumn[3] = new std::string("$\\sigma_{0}$");

    // Filas: cada fila representa un valor de Ex
    for (int i = 0; i < 4; i++) {
        file << *firstColumn[i] << ",";
        for (int j = 0; j < 2; j++) {
            file <<  TString::Format("$\\num{%.12f(%.10f)}$",sigma[j][i],usigma[j][i]);
            if (j < 4) file << ","; // no agregar coma al final de línea
        }
        file << "\n";
    }

    file.close();
    std::cout << "Archivo sigmas.csv exportado correctamente.\n";
}

TH1* read_tree(double Ex, TCanvas* c, int squareCanvas,THStack* hs,double& sigmaOut, double& usigmaOut,EColor color = kBlue, int incIdx=0) {

    ROOT::EnableImplicitMT(); // Enable multithreading
    ROOT::RDataFrame df {"SimulationTree", TString::Format("./Outputs/tree_Ex_%.2f_incIdx%i.root", Ex,incIdx)};
    df.Describe().Print();
    
    auto hEx {df.Histo1D({"hEx", TString::Format("Rec Ex=%.2f;E_{x} [MeV];Counts", Ex), 200, -2.5, 3}, "exrec")};
    
    // Clona el histograma
    auto* hClone = (TH1*) hEx->Clone();

    // Asigna color de relleno (aquí ejemplo simple)
    static int colorIndex =2;  // Empieza en 2 (rojo)
    //hClone->SetFillColor(colorIndex++);

    hClone->SetLineColor(color);  // Borde negro opciona
    hClone->SetLineWidth(3);
    // Añade al THStack
    hs->Add(hClone);  

    // Dibujamos el gráfico por separado
    c->cd(squareCanvas);

    // Ajustamos a una Voigt -> Convolución gauss+ briet-wigner
    double xmin=Ex-5;
    double xmax=Ex+5;

    TF1 *fitVoigt = new TF1("voigt", [](double *x, double *par) {
        return par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]);
    }, xmin,xmax, 4);

    // Introducimos estimaciones iniciales para obtener las Voigt ajustadas: 
    double sigmaInicial{0.05};
    double gammaInicial{0.1};
    if (incIdx == 0 || incIdx == 2){
        sigmaInicial=0.2;
    }

    fitVoigt->SetParameters(10e6, Ex, sigmaInicial, gammaInicial); // estimación inicial
    fitVoigt->SetLineColor(color-5); 
    fitVoigt->SetNpx(5000); // o más, por defecto es 100
    fitVoigt->SetLineWidth(4);
    hClone->Fit(fitVoigt, "RM");
    

    
    hClone->Draw();
    fitVoigt->Draw("same");

    
    
    // Parámetros:
    // par[0] = Amplitud
    // par[1] = Media (μ)
    // par[2] = Sigma (anchura gaussiana)
    // par[3] = Gamma (anchura Lorentziana)
    double sigma = fitVoigt->GetParameter(2);
    double usigma = fitVoigt->GetParError(2);

    sigmaOut=sigma;
    usigmaOut=usigma;

    return hClone;
}
// Dibuja sin apilar (nostack)

void Read_tree_aux(int incIdx, double sigma[2],double usigma[2]){

    gStyle->SetOptStat(0);  // Desactiva globalmente
    auto* c0 {new TCanvas {"c0", "c0"}};
    auto* c1 {new TCanvas {"c1", "c1"}};
    auto* hs = new THStack("hs", ";E_{x} [MeV];Counts");
    c0->DivideSquare(4);

    // Guardamos los histogramas 


    std::vector<std::pair<TH1*, double>> hList;

    hList.push_back({read_tree(0.0, c0, 1, hs, sigma[0],usigma[0] ,kRed, incIdx), 0.0});
    hList.push_back({read_tree(0.20, c0, 2, hs,sigma[1],usigma[1] ,kBlue, incIdx), 0.20});
    //hList.push_back({read_tree(0.47, c0, 3, hs,sigma[2] ,kGreen,  incIdx), 0.47});


    // Crea la leyenda
    auto* legend = new TLegend(0.7, 0.7, 0.9, 0.9);  // (xmin, ymin, xmax, ymax)
    auto* hSum {new TH1D {"hEx", ";E_{x} [MeV];Counts", 200, -2.5, 3}};

    for (auto& [h, ex] : hList) {
        legend->AddEntry(h, TString::Format("E_{x}=%.2f MeV", ex), "f");  // "f" = fill 
        hSum->Add(h);
    }
    

    c1->cd();
    c1->DivideSquare(2);

    c1->cd(1);
    
    hs->Draw("nostack");
    hs->GetXaxis()->SetRangeUser(-0.3, 0.6);  // Solo muestra de -0.2 a 0.6 en X
    legend->Draw();


    // Antes de modificar
    std::vector<TH1*> cleanHistos;

    TIter next(hs->GetHists());
    TH1* h = nullptr;

    while ((h = (TH1*)next())) {
        TH1* hclone = (TH1*) h->Clone();  // Crea copia
        hclone->GetListOfFunctions()->Clear();  // Borra ajustes
        cleanHistos.push_back(hclone);
    }
    THStack* hs_clean = new THStack("hs_clean", "");  // no importa el título aquí
    for (auto& h : cleanHistos) hs_clean->Add(h);
    // Dibujamos el HStack

    c1->cd(2);

    hSum->SetLineColor(kGreen-5);
    hSum->SetLineWidth(3);


    while ((h = (TH1*)next())) {
        h->GetListOfFunctions()->Clear();  // Elimina todos los TF1 asociados
    }


    hSum->Draw("hist");  // sin "same" porque es primero
    hSum->GetXaxis()->SetRangeUser(-0.3, 0.6);  // Solo muestra de -0.2 a 0.6 en X
    hs_clean->Draw("nostack same");

    legend->AddEntry(hSum, "Apilada", "f");
    legend->Draw();

    // Guardamos
    c0->Update();
    
    //c0->SaveAs(TString::Format("Graficas/ExHisto/Rec_incIdx%i.pdf",incIdx));
    c1->SaveAs(TString::Format("Graficas/ExHisto/Rec_incIdx%i_single.pdf",incIdx));

}


void Read_tree(){

    int uncertaintySelector[4] = {0,1,2,3};
    double sigma[2][4] {0};
    double sigmaIdx[2] {0};
    double usigma[2][4] {0};
    double usigmaIdx[2] {0};

    for (int i = 0; i < 4; i++) {
        Read_tree_aux(uncertaintySelector[i],sigmaIdx,usigmaIdx);
        for (int j=0; j<2; j++){
            sigma[j][i]=sigmaIdx[j];
            usigma[j][i]=usigmaIdx[j];
        }
    }
    ExportSigmasToCSV(sigma,usigma,"SigmasTab.csv");

}