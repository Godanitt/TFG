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
void xs()
{
    auto *gs{new TGraphErrors{"Inputs/xs/s12_p1i.dat", "%lg %lg"}};
    gs->SetTitle("s_{1/2}");
    auto *gp{new TGraphErrors{"Inputs/xs/p12_p1i.dat", "%lg %lg"}};
    gp->SetTitle("p_{1/2}");

    auto *mg{new TMultiGraph};
    mg->SetTitle(";#theta_{CM} [#circ];d#sigma/d#Omega [mb/sr]");

    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9); // esquina superior derecha

    std::vector<std::string> name = {"Ex=0, s12", "Ex=0.2, p12"};
    int i{1};

    std::vector<TGraphErrors *> graphs = {gs, gp};
    for (size_t j = 0; j < graphs.size(); ++j)
    {
        graphs[j]->SetLineWidth(2);
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

    auto *c0{new TCanvas{"c0", "xs canvas"}};
    // c0->DivideSquare(2);
    // c0->cd(1);

    mg->Draw("apl");
    leg->Draw();
    // c0->cd(2);
    // hThetaCM->Draw();

    c0->SaveAs(TString::Format("Graficas/Seccion_Eficaz.pdf"));
}
