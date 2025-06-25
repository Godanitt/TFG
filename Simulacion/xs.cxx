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
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"

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
}
