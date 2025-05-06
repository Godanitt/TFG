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
double ApplyAngleResolution(double theta4old, double FWHM)
{
    double sigmatheta{FWHM  / 2.55};
    return gRandom->Gaus(theta4old, sigmatheta);
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
                                  //verticeCentrado.Z() * 10 + 175
                                  135
                                };

    return verticeUsual;
}