#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"
#include "TMath.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "TVector.h"


struct ResultadoTriggerL1
{
    double range;
    bool hasStoppedInTPC;      

};

ResultadoTriggerL1 TriggerL1(double t3, double theta3, ROOT::Math::XYZPoint vertex, ROOT::Math::XYZVector direction, ActPhysics::SRIM &srim) {

    // Actar
    ActRoot::TPCParameters tpc {"Actar"};
    // Tamanhos en mm:
    auto xActar {tpc.X()}; // o mesmo para Y,Z
    auto yActar {tpc.X()}; // o mesmo para Y,Z
    auto zActar {tpc.X()}; // o mesmo para Y,Z
    
    bool hasStoppedinTPC {false};

    double rango{srim.EvalRange("light", t3)};
    ROOT::Math::XYZPoint projectedRange {vertex + rango*direction};

    double x = projectedRange.X();
    double y = projectedRange.Y();
    double z = projectedRange.Z();

    if (x<xActar && y<yActar && z<zActar){

        hasStoppedinTPC = true;
    }
    ResultadoTriggerL1 resTriggerL1;

    resTriggerL1.range=rango;
    resTriggerL1.hasStoppedInTPC=hasStoppedinTPC;

    return resTriggerL1;    
}