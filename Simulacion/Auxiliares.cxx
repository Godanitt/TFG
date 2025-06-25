#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"

#include "TVector3.h" 
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
//#######################################################
/*
 Aplicamos resolucion a la variable thetalab. Entrada:
 - double thetalab, 

 Salida: 
 - double thetalab' (con resolucion)
*/
double ApplyAngleResolution(double thetaold, double FWHM)
{
    double sigmatheta{FWHM  / 2.55};
    return gRandom->Gaus(thetaold, sigmatheta);
}


//#######################################################
/*
 Definimos el vértice centrado en el detector  y con coordenadas en mm. Entrada: 
 - double xACTAR, 
 - double yACTAR, 
 - double zACTAR 
*/
ROOT::Math::XYZPoint SampleVertex(double xACTAR, double yACTAR, double zACTAR)
{
    ROOT::Math::XYZPoint verticeCentrado{
    gRandom->Uniform(-xACTAR / 2, xACTAR / 2) * 0.1,
    gRandom->Gaus(0, 5) * 0.1,
    gRandom->Gaus(0, 5) * 0.1};

    ROOT::Math::XYZPoint verticeUsual{verticeCentrado.X() * 10 + xACTAR / 2,
                                  verticeCentrado.Y() * 10 + yACTAR / 2,
                                  verticeCentrado.Z() * 10 +  zACTAR / 2
                                };

    return verticeUsual;
}

//####################################################### 
/*
La entrada es i. Si i=1 -> Se elimina el título.
*/
void SetMyStyle() {
    // Elimina el panel de estadísticas
    gStyle->SetOptStat(0);

    // Margenes: 
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.2);
    gStyle->SetPadTopMargin(0.08);

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



}