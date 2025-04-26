#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"

TH1* read_tree(double Ex, TCanvas* c, int squareCanvas,THStack* hs,EColor color = kBlue) {

    ROOT::EnableImplicitMT(); // Enable multithreading
    ROOT::RDataFrame df {"SimulationTree", TString::Format("./Outputs/tree_Ex_%.2f.root", Ex)};
    df.Describe().Print();
    
    auto hEx {df.Histo1D({"hEx", TString::Format("Rec Ex=%.2f;E_{x} [MeV];Counts", Ex), 200, -3, 7}, "exrec")};
    
    // Clona el histograma
    auto* hClone = (TH1*) hEx->Clone();

    // Asigna color de relleno (aquí ejemplo simple)
    static int colorIndex = 2;  // Empieza en 2 (rojo)
    hClone->SetFillColor(colorIndex++);
    hClone->SetLineColor(kBlack);  // Borde negro opcional
    
    // Añade al THStack
    hs->Add(hClone);  

    // Dibujamos el gráfico por separado
    c->cd(squareCanvas);
    hClone->Draw();
    return hClone;
}
// Dibuja sin apilar (nostack)


void Read_tree(){
    auto* c0 {new TCanvas {"c0", "c0"}};
    auto* hs = new THStack("hs", "THStack Eex;Ex [MeV];Cuentas");
    c0->DivideSquare(4);

    // Guardamos los histogramas 

    std::vector<std::pair<TH1*, double>> hList;

    hList.push_back({read_tree(0.0, c0, 1, hs,kRed), 0.0,});
    hList.push_back({read_tree(0.2, c0, 2, hs,kBlue), 0.2,});
    hList.push_back({read_tree(0.4, c0, 3, hs,kGreen), 0.4,});
  

    // Dibujamos el HStack
    c0->cd(4);
    hs->Draw("nostack");

    // Crea la leyenda
    auto* legend = new TLegend(0.7, 0.7, 0.9, 0.9);  // (xmin, ymin, xmax, ymax)
    for (auto& [h, ex] : hList) {
        legend->AddEntry(h, TString::Format("E_{x}=%.2f MeV", ex), "f");  // "f" = fill style
    }

    legend->Draw();
    // Guardamos
    c0->Update();
    c0->SaveAs("Graficas/ExRec.pdf");
}