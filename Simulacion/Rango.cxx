#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "ActCrossSection.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include "TCanvas.h"

#include "TVector3.h" 
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "TLine.h"



/*
Aquí vamos a relacionar rango-energía-angulo, obteniendo las gráficas 
*/
TGraph* Rango(ActPhysics::Kinematics *kin, ActPhysics::SRIM &srim) {

auto* g3 {kin->GetKinematicLine3()};

int n = g3->GetN();

std::vector<double> thetas, rangos,rangosReally, energias;

double rango;
double rangoReally;

for (int i=0; i<n; ++i) {
    double theta, energia;
    g3->GetPoint(i, theta, energia);


    // Si usas la función fRange(E):
    rango =(2.0/4.0) *  srim.EvalRange("lightInSil", energia);

    rangoReally = (1.0/1.0) * srim.EvalRange("lightInSil", energia);
    // Se multiplica por 1/2 porque no sabemos cuanto afecta que el rango sea en XY
    // hay que ponderar ya que range=sqrt(x**2+y**2+z**2) -> si rXY**2=range**2-z**2



    // O si usas un TGraph gRangeE con rango(E):
    // rango = gRangeE->Eval(energia);

    thetas.push_back(theta);
    rangos.push_back(rango);
    rangosReally.push_back(rangoReally);
    energias.push_back(energia);
}

// Crear el nuevo TGraph rango vs theta
TGraph* gRangoVsTheta = new TGraph(n, &thetas[0], &rangos[0]);
TGraph* gRangoVsE = new TGraph(n, &energias[0], &rangosReally[0]);




auto *c10{new TCanvas{"c10", "c10", 440,440}};
gRangoVsTheta->SetLineWidth(3);
gRangoVsTheta->GetXaxis()->SetTitle("#theta [ ^{#circ} ]");  // eje X
gRangoVsTheta->GetYaxis()->SetTitle("Rango [mm]");        // eje Y
gRangoVsTheta->SetTitle("");        // eje Y
gRangoVsTheta->Draw("al");
gRangoVsTheta->GetYaxis()->SetRangeUser(0, 200);

c10->SetRightMargin(0.05);
c10->SetLeftMargin(0.16);
c10->SetBottomMargin(0.15);
c10->SaveAs(TString::Format("/home/daniel/GitHub/TFG/Memoria/Imagenes/RangoThetaTeo.pdf"));



auto *c11{new TCanvas{"c11", "c11", 440,440}};
gRangoVsE->SetLineWidth(3);
gRangoVsE->GetXaxis()->SetTitle("E [MeV]"); // eje X
gRangoVsE->GetYaxis()->SetTitle("Rango [mm]");        // eje Y
gRangoVsE->SetTitle("");        // eje Y
gRangoVsE->Draw("al");
c11->SetRightMargin(0.05);
c11->SetLeftMargin(0.16);
c11->SetBottomMargin(0.15);
c11->SaveAs(TString::Format("/home/daniel/GitHub/TFG/Memoria/Imagenes/RangoEnergiaTeo.pdf"));




return gRangoVsTheta;
}