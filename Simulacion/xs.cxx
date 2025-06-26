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


void ExportSigmasToCSV(const double sigma[2][2], const std::string& filename = "SeccionEficaz.csv") {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "No se pudo abrir el archivo para escritura.\n";
        return;
    }

    // Cabecera (opcional)
    file << ",$\\sigma(0.0)$ [mb] , $\\sigma(0.20)$ [mb] \n";

    std::string* firstColumn[4];  // array de 4 punteros a std::string

    firstColumn[0] = new std::string("Trapecio");
    firstColumn[1] = new std::string("Simpson");

    // Filas: cada fila representa un valor de Ex
    for (int i = 0; i < 2; i++) {
        file << *firstColumn[i] << ",";
        for (int j = 0; j < 2; j++) {
            file <<  TString::Format("$\\num{%.4f}$",sigma[j][i]);
            if (j < 2) file << ","; // no agregar coma al final de línea
        }
        file << "\n";
    }
}

void ExportsCSVToTex(const std::string& filenameCSV,const std::string& filenameTEX) {

    std::ifstream csv(filenameCSV);
    std::ofstream tex("/home/daniel/GitHub/TFG/Memoria/Cuerpo/"+filenameTEX);

    tex << "\\begin{tabular}{llllllll} \\hline\n";  // ajusta columnas según el CSV

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

}

void xs()
{
    auto *gs{new TGraphErrors{"Inputs/xs/s12_p1i.dat", "%lg %lg"}};
    gs->SetTitle("s_{1/2}");
    auto *gp{new TGraphErrors{"Inputs/xs/p12_p1i.dat", "%lg %lg"}};
    gp->SetTitle("p_{1/2}");

    auto *mg{new TMultiGraph};
    mg->SetTitle(";#theta_{CM} [#circ];d#sigma/d#Omega [mb/sr]");

    TLegend *leg = new TLegend(0.6, 0.72, 0.99, 0.99); // esquina superior derecha

    std::vector<std::string> name = {"E_{ex}=0.0 MeV, 2s_{1/2}", "E_{ex}=0.2 MeV, 1p_{1/2}"};
    int i{1};

    std::vector<TGraphErrors *> graphs = {gs, gp};
    for (size_t j = 0; j < graphs.size(); ++j)
    {
        graphs[j]->SetLineWidth(4);
        graphs[j]->SetLineColor(j + 1); // 1 = negro, 2 = rojo
        leg->AddEntry(graphs[j], name[j].c_str(), "l");
        mg->Add(graphs[j]);
    }

    // Sampling
    ActSim::CrossSection xs{};
    xs.ReadFile("Inputs/xs/p12_p1i.dat");
    xs.Draw();

    auto *hThetaCM{new TH1D{"hThetaCM", "CM;#theta_{CM};Counts", 300, 0, 180}};
    for (int i = 0; i < 10000000; i++)
    {
        hThetaCM->Fill(xs.SampleHist());
    }


    // Margenes: 
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.01);
    gStyle->SetPadTopMargin(0.01);
    gStyle->SetPadBottomMargin(0.12);

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
    leg->SetTextFont(62);     // Fuente en negrita
    auto *c0{new TCanvas{"c0", "xs canvas"}};
    // c0->DivideSquare(2);
    // c0->cd(1);

    mg->Draw("apl");
    leg->Draw();
    // c0->cd(2);
    // hThetaCM->Draw();
    

    c0->SaveAs(TString::Format("/home/daniel/GitHub/TFG/Memoria/Imagenes/Seccion_Eficaz.pdf"));

    // Ahora vamos a integrar 
    // ========== INTEGRACIÓN MÉTODO 1: NUMÉRICO ==========

    double sigma[2][2] {0};

    for (size_t j = 0; j < graphs.size(); ++j)
    {
        double integral1 = 0;
        int N = graphs[j]->GetN();
        double *x = graphs[j]->GetX();
        double *y = graphs[j]->GetY();

        for (int i = 0; i < N - 1; ++i)
        {
            double theta1 = x[i] * TMath::DegToRad();
            double theta2 = x[i + 1] * TMath::DegToRad();
            double dtheta = theta2 - theta1;

            // Promedio de valores
            double dsdo1 = y[i];
            double dsdo2 = y[i + 1];
            double avg = 0.5 * (dsdo1 * sin(theta1) + dsdo2 * sin(theta2));

            integral1 += avg * dtheta;
        }

        integral1 *= 2 * TMath::Pi(); // integrar en dΩ = 2π sinθ dθ
        std::cout << "[Método 1] Integral total (" << name[j] << "): " << integral1 << " mb" << std::endl;
        sigma[j][0]=integral1;
    }
    for (size_t j = 0; j < graphs.size(); ++j)
    {
    double integralSimpson = 0;
    int N = graphs[j]->GetN();

    // Asegurarse de que haya un número impar de puntos (par de intervalos)
    if ((N - 1) % 2 != 0)
    {
        std::cout << "Advertencia: el número de puntos no es válido para Simpson, se omite el último punto.\n";
        N -= 1;
    }

    double *x = graphs[j]->GetX();
    double *y = graphs[j]->GetY();

    for (int i = 0; i < N - 1; i += 2)
    {
        double theta0 = x[i] * TMath::DegToRad();
        double theta1 = x[i + 1] * TMath::DegToRad();
        double theta2 = x[i + 2] * TMath::DegToRad();

        double h = (theta2 - theta0) / 2.0;

        double f0 = y[i] * sin(theta0);
        double f1 = y[i + 1] * sin(theta1);
        double f2 = y[i + 2] * sin(theta2);

        integralSimpson += (h / 3.0) * (f0 + 4.0 * f1 + f2);
    }

    integralSimpson *= 2.0 * TMath::Pi();  // para cubrir dΩ = 2π sinθ dθ

    std::cout << "[Método 3: Simpson] Integral total (" << name[j] << "): "
              << integralSimpson << " mb" << std::endl;

    sigma[j][1]=integralSimpson;

    }

    ExportSigmasToCSV(sigma,"SeccionEficaz.csv");

    ExportsCSVToTex("SeccionEficaz.csv","SeccionEficaz.tex");

    double alpha1{6000.0*6.0*24.0*3600.0*0.7826*4.9945*((1e-8))*25.5*(sigma[0][0])};
    double alpha2{6000.0*6.0*24.0*3600.0*0.7826*4.9945*((1e-8))*25.5*(sigma[1][0])};

    std::ofstream tex("/home/daniel/GitHub/TFG/Memoria/Cuerpo/alphas.tex");
    std::cout<<"sigma="<<sigma[0][0]<<"\n";
    tex<<"\\alpha (0.0)= \\frac{"<<alpha1<<"}{N_{\\text{inter}}} \\qquad  \n";
    tex<<"\\alpha (0.2)=\\frac{"<<alpha2<<"}{N_{\\text{inter}}}  \n";

    }
