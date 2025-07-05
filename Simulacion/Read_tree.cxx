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

void ExportSigmasToCSV(const double sigma[2][5],const double usigma[2][5], const std::string& filename = "sigmas.csv", const std::string& nombre = "sigma") {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "No se pudo abrir el archivo para escritura.\n";
        return;
    }

    // Cabecera (opcional)
    file << ",$\\" + nombre +"(0.0)$ [keV] , $\\" + nombre+ "(0.20)$ [keV] \n"; 

    std::string* firstColumn[4];  // array de 4 punteros a std::string

    firstColumn[0] = new std::string("$\\"+nombre+"_{tot}$");
    firstColumn[1] = new std::string("$\\"+nombre+"_{str}$");
    firstColumn[2] = new std::string("$\\"+nombre+"_{\\theta}$");
    firstColumn[3] = new std::string("$\\"+nombre+"_{0}$");
    // Filas: cada fila representa un valor de Ex
    for (int i = 0; i < 4; i++) {
        file << *firstColumn[i] << ",";
        for (int j = 0; j < 2; j++) {
            file <<  TString::Format("$\\num{%.1f(%.1f)}$",1000.0*sigma[j][i],1000.0*usigma[j][i]);
            if (j < 4) file << ","; // no agregar coma al final de línea
        }
        file << "\n";
    }
}

TH1* read_tree(double Ex, TCanvas* c, int squareCanvas,THStack* hs,double& sigmaOut, double& usigmaOut,double& gammaOut, double& ugammaOut, EColor color = kBlue, int incIdx=0) {

    ROOT::EnableImplicitMT(); // Enable multithreading
    ROOT::RDataFrame df {"SimulationTree", TString::Format("./Outputs/tree_Ex_%.2f_incIdx%i.root", Ex,incIdx)};
    df.Describe().Print();
    
    auto hEx {df.Histo1D({"hEx", TString::Format("Rec Ex=%.2f;E_{x} [MeV];Cuentas", Ex), 1200, -10, 10}, "exrec")};
    
   // if (incIdx==4) {   
   //     hEx=df.Histo1D({"hEx", TString::Format("Rec Ex=%.2f;E_{x} [MeV];Cuentas", Ex), 6000, -0.2 ,0.4}, "exrec");
   // }


    double factorHisto{0};

    double interacciones{5000000};

    if (Ex==0.0) {factorHisto=336920/interacciones;}

    if (Ex==0.2) {factorHisto=259406/interacciones;}

    
    hEx->Scale(factorHisto);   


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
    double xmin=Ex-10;
    double xmax=Ex+10;


    double sigma{0};
    double usigma{0};
    double gamma{0};
    double ugamma{0};

    if (incIdx!=4){
        TF1 *fitVoigt = new TF1("voigt", [](double *x, double *par) {
            return par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]);
        }, xmin,xmax, 4);
        // Introducimos estimaciones iniciales para obtener las Voigt ajustadas: 
        double sigmaInicial{0.00};
        double gammaInicial{0.3};
        if (Ex==0.0)
            double gammaInicial=0.1;
        if (incIdx == 0 || incIdx == 2){
            sigmaInicial=0.2;
        }


        fitVoigt->SetParameters(10e4, Ex, sigmaInicial, gammaInicial); // estimación inicial
        fitVoigt->SetLineColor(color-5); 
        fitVoigt->SetNpx(5000); // o más, por defecto es 100
        fitVoigt->SetLineWidth(4);
        hClone->Fit(fitVoigt, "RM");
        

        
        hClone->Draw("hist E0");
        fitVoigt->Draw("same");

        
        
        // Parámetros:
        // par[0] = Amplitud
        // par[1] = Media (μ)
        // par[2] = Sigma (anchura gaussiana)
        // par[3] = Gamma (anchura Lorentziana)
        sigma = fitVoigt->GetParameter(2);
        usigma = fitVoigt->GetParError(2);
        gamma = fitVoigt->GetParameter(3);
        ugamma = fitVoigt->GetParError(3);
    }

    if (incIdx==4){
        TF1 *fitGauss = new TF1("gauss", [](double *x, double *par) {
            return par[0] * TMath::Gaus(x[0], par[1], par[2], true); // true = normalizada
        }, xmin, xmax, 3);

        // Estimaciones iniciales de parámetros:
        double sigmaInicial{0.05};
        if (incIdx == 0 || incIdx == 2){
            sigmaInicial = 0.2;
        }



        // Ajuste inicial: par[0]=Amplitud, par[1]=Media, par[2]=Sigma
        fitGauss->SetParameters(10e6, Ex, sigmaInicial); 
        fitGauss->SetLineColor(color - 5);
        fitGauss->SetNpx(5000); 
        fitGauss->SetLineWidth(4);

        // Ajuste al histograma
        hClone->Fit(fitGauss, "RM"); 

        // Dibujo
        hClone->Draw("hist E0");
        fitGauss->Draw("same");

        // Parámetros:
        // par[0] = Amplitud
        // par[1] = Media (μ)
        // par[2] = Sigma (anchura Gaussiana)
        sigma = fitGauss->GetParameter(2);
        usigma = fitGauss->GetParError(2);

    }

    sigmaOut=sigma;
    usigmaOut=usigma;
    gammaOut=gamma;
    ugammaOut=ugamma;



    return hClone;
}
// Dibuja sin apilar (nostack)

void Read_tree_aux(int incIdx, double sigma[2],double usigma[2],double gamma[2],double ugamma[2]){

    gStyle->SetOptStat(0);  // Desactiva globalmente
    auto* c0 {new TCanvas {"c0", "c0"}};
    auto* c1 {new TCanvas {"c1", "c1",750, 300}};
    auto* hs = new THStack("hs", ";E_{x} [MeV];Cuentas");
    c0->DivideSquare(4);

    // Guardamos los histogramas 


    std::vector<std::pair<TH1*, double>> hList;

    hList.push_back({read_tree(0.0, c0, 1, hs, sigma[0],usigma[0],gamma[0],ugamma[0],kRed, incIdx), 0.0});
    hList.push_back({read_tree(0.20, c0, 2, hs,sigma[1],usigma[1],gamma[1],ugamma[1],kBlue, incIdx), 0.20});
    //hList.push_back({read_tree(0.47, c0, 3, hs,sigma[2] ,kGreen,  incIdx), 0.47});


    // Crea la leyenda
    auto* legend = new TLegend(0.64, 0.75, 0.95, 0.95);  // (xmin, ymin, xmax, ymax)
    auto* hSum {new TH1D {"hEx", ";E_{x} [MeV];Cuentas", 1200, -10, 10}};

    for (auto& [h, ex] : hList) {
        legend->AddEntry(h, TString::Format("E_{x}=%.2f MeV", ex), "f");  // "f" = fill 
        hSum->Add(h);
    }
    
    gStyle->SetLabelSize(0.05,"XYZ");  // Tamaño de los números en el eje X
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetTitleSize(0.05,"XYZ");  // Tamaño de los números en el eje X
    gStyle->SetTitleFont(62);      // Sin segundo argumento → título general
    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetTitleSize(0.05); 
    // Margenes: 
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.14);

    // Grosor de los ejes
    gStyle->SetLineWidth(2);      // Más grueso el eje X
    gStyle->SetLineWidth(2);      // Más grueso el eje Y


    c1->cd();
    c1->DivideSquare(2);

    c1->cd(1);
    
    hs->Draw("lh nostack");
    hs->GetXaxis()->SetRangeUser(-0.4, 0.8);  // Solo muestra de -0.2 a 0.6 en X
    
    legend->SetTextSize(0.05);  // Tamaño del texto
    legend->SetTextFont(62);     // Fuente en negrita
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
    hs_clean->Draw("hist nostack same");

    legend->AddEntry(hSum, "Suma", "f");
    legend->Draw();

    // Guardamos
    c0->Update();
    
    //c0->SaveAs(TString::Format("Graficas/ExHisto/Rec_incIdx%i.pdf",incIdx));
    c1->SaveAs(TString::Format("/home/daniel/GitHub/TFG/Memoria/Imagenes/Rec_incIdx%i_single.pdf",incIdx));

}


void Read_tree(){

    int uncertaintySelector[5] = {0,1,2,3,4};
    double sigma[2][5] {0};
    double sigmaIdx[2] {0};
    double usigma[2][5] {0};
    double usigmaIdx[2] {0};

    double gamma[2][5] {0};
    double gammaIdx[2] {0};
    double ugamma[2][5] {0};
    double ugammaIdx[2] {0};

    for (int i = 0; i < 4; i++) {
        Read_tree_aux(uncertaintySelector[i],sigmaIdx,usigmaIdx,gammaIdx,ugammaIdx);
        for (int j=0; j<2; j++){
            sigma[j][i]=sigmaIdx[j];
            usigma[j][i]=usigmaIdx[j];
            gamma[j][i] = gammaIdx[j];       // CORREGIDO
            ugamma[j][i] = ugammaIdx[j];     // CORREGIDO
            std::cout<<"sigma"<<sigmaIdx[j]<<"\n";
            std::cout<<"gamma"<<gammaIdx[j]<<"\n";
        }
    }
    ExportSigmasToCSV(sigma,usigma,"SigmasTab.csv");
    ExportSigmasToCSV(gamma,ugamma,"GammaTab.csv","Gamma");



    std::cout << "Archivo sigmas.csv exportado correctamente.\n";

    std::ifstream csv("SigmasTab.csv");
    std::ofstream tex("/home/daniel/GitHub/TFG/Memoria/Cuerpo/SigmasTab.tex");

    tex << "\\begin{tabular}{llll} \\hline\n";  // ajusta columnas según el CSV

    tex << "\\toprule \n";  // ajusta columnas según el CSV

    std::string linea;
    bool primera = true;
    bool isPrimero=true;
    while (std::getline(csv, linea)) {
        std::istringstream ss(linea);
        std::string valor;
        bool primero = true;

        

        while (std::getline(ss, valor, ',')) {
            if (!primero) tex << " & ";
            tex << valor;
            primero = false;
        }
        if (isPrimero==true){tex << " \\\\ \\midrule \n";}
        else {tex << "\\\\ \n";}

        isPrimero=false;

    }
    tex << " \\bottomrule \n";
    tex << "\\end{tabular}\n";

    std::cout << "Archivo gammas.csv exportado correctamente.\n";

    isPrimero=true;

    std::ifstream csv2("GammaTab.csv");
    std::ofstream tex2("/home/daniel/GitHub/TFG/Memoria/Cuerpo/GammaTab.tex");

    tex2 << "\\begin{tabular}{llll} \\hline\n";  // ajusta columnas según el CSV

    tex2 << "\\toprule \n";  // ajusta columnas según el CSV

    while (std::getline(csv2, linea)) {
        std::istringstream ss(linea);
        std::string valor;
        bool primero = true;

        while (std::getline(ss, valor, ',')) {
            if (!primero) tex2 << " & ";
            tex2 << valor;
            primero = false;
        }
        if (isPrimero==true){tex2 << " \\\\ \\midrule \n";}
        else {tex2 << "\\\\ \n";}

        isPrimero=false;

    }
    tex2 << " \\bottomrule \n";


    tex2 << "\\end{tabular}\n";

    
}