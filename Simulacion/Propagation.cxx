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

// #######################################################
/*
 Aplicamos resolucion relacionada con las pérdidas de enerǵia en el gas. Entrada:
 - SRIM::*srim (puntero) -> Nos dincIdxice el gas en el que está y como van actuar las pérdidas
 - which: tipo de gas (beam, light, heavy)
 - double Tini: energía inicial (en MeV)
 - double dist: distancia recorrida (en mm)

 Salida:
 - double tNew' (con resolucion)
*/
double ApplyStraggling(ActPhysics::SRIM *srim, std::string which, double tIni, double dist)
{
    double Rini{srim->EvalRange(which, tIni)};
    double uini{srim->EvalLongStraggling(which, Rini)};

    double Rleft{Rini - dist};
    if (Rleft < 0)
        return -1;
    double uleft{srim->EvalLongStraggling(which, Rleft)};

    double udist{TMath::Sqrt(TMath::Power(uini, 2) - TMath::Power(uleft, 2))};

    double d{gRandom->Gaus(dist, udist)};
    double Rleftnew{Rini - d};

    double tNew{srim->EvalEnergy(which, Rleftnew)};

    return tNew;
}
// #######################################################
/*
 Aplicamos resolucion a la variable energía. Entrada:
 - double deltaE,

 Salida:
 - double deltaE' (con resolucion)
*/
double ApplySilResolution(double deltaE)
{
    return gRandom->Gaus(deltaE, (0.0213 * TMath::Sqrt(deltaE)) / (2.35));
}
// #######################################################

struct Resultado
{
    // Perdidas de energia
    double dt3sil0;
    double dt3sil1;

    // Distancias en las que se perdieron la energía:
    double sildist0;
    double sildist1;

    // Puntos de imapacto del silicio:silIndex1
    ROOT::Math::XYZPoint silPoint1;

    // Booleano que nos dice si se ha medido en el layer 0:
    bool isInSil0;

    // Booleano que nos dice si se ha medido en el layer 1:
    bool isInSil1;

    // Booleano que nos dice si se ha paratdo antes del silicio 0
    bool stoppedBeforeSil0;

    // Valor del rango que alcanza la partícula antes de pararase
    double rangeBeforeSil0;
   

};

Resultado Propagation(double t3, double theta3, ROOT::Math::XYZPoint vertex, ROOT::Math::XYZVector direction, ActPhysics::SRIM &srim, ActPhysics::SilSpecs *sils, int incIdx)
{
    Resultado res;
      int silIndex0 {};
        ROOT::Math::XYZPoint silPoint0 {};
        std::string layer0 {};
        for(const auto& layer : {"f0", "l0", "r0"})
        // for(const auto& layer : {"f0"})
        {
            auto [idx, point] = sils->FindSPInLayer(layer, vertex, direction);
            if(idx != -1)
            {
                silIndex0 = idx;
                silPoint0 = point;
                layer0 = layer;
                break;
            }
        }
        // Continue if no silicon is reached geometrically

        // Propagate to it
        auto dist0 {(vertex - silPoint0).R()};
        auto TAtSil0 {srim.Slow("light", t3, dist0)};

        if (incIdx==0 ||  incIdx==1 ){        
            TAtSil0=srim.SlowWithStraggling("light", t3, dist0);
        }

        ///////////// SILICON 0 ///////////////////
        // Angle with normal
        auto thick0 {sils->GetLayer(layer0).GetUnit().GetThickness()};
        auto normal {sils->GetLayer(layer0).GetNormal()};
        auto angleNormal0 {TMath::ACos(direction.Unit().Dot(normal.Unit()))};
        // Eloss
        auto TAfterSil0 {srim.Slow("lightInSil", TAtSil0, thick0, angleNormal0)};
        auto eloss0 {TAtSil0 - TAfterSil0};

        if (incIdx==0 ||  incIdx==1 ){
            TAfterSil0=srim.SlowWithStraggling("lightInSil", TAtSil0, thick0, angleNormal0);
            eloss0 = ApplySilResolution(TAtSil0 - TAfterSil0) ;
        }

        // std::cout << "TAtSil0 : " << TAtSil0 << '\n';
        // std::cout << "thick0 : " << thick0 << '\n';
        // std::cout << "angleNormal0 : " << angleNormal0 << '\n';
        // std::cout << "eloss0 : " << eloss0 << '\n';

        //////////// SILICON 1 ////////////////////
        // only if layer0 == "f0"
        int silIndex1 {};
        ROOT::Math::XYZPoint silPoint1 {};
        double dist1 {};
        double TAtSil1 {};
        double eloss1 {};
        double TAfterSil1 {};
        if(layer0 == "f0" && TAfterSil0 > 0)
        {
            // Propagate to sil1
            std::tie(silIndex1, silPoint1) = sils->FindSPInLayer("f1", vertex, direction);
            // Eloss in intergas = region between f0 and f1
            dist1 = (silPoint0 - silPoint1).R();
            TAtSil1 = srim.Slow("light", TAfterSil0, dist1);
            
            if (incIdx==0 ||  incIdx==1 ){        
                TAtSil1=srim.SlowWithStraggling("light", TAfterSil0, dist1);
            }

            auto thick1 {sils->GetLayer("f1").GetUnit().GetThickness()};
            auto normal1 {sils->GetLayer("f1").GetNormal()};
            auto angleNormal1 {TMath::ACos(direction.Unit().Dot(normal.Unit()))};
            TAfterSil1 = srim.Slow("lightInSil", TAtSil1, thick1, angleNormal1);
            eloss1 = TAtSil1 - TAfterSil1;

            if (incIdx==0 ||  incIdx==1  ){
                TAfterSil1=srim.SlowWithStraggling("lightInSil", TAtSil1, thick1, angleNormal1);
                eloss1 = ApplySilResolution(TAtSil1 - TAfterSil1) ;
            }
        }

    //%%%%%%%%%%%layer0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Guardamos los datos: 


    bool isInSil0{silIndex0 != -1 && 0 < TAtSil0 && !(TAfterSil0>0)};
    bool isInSil1{!(TAtSil1 <= 0) && silIndex0 != -1  && silIndex1!=-1 &&  !(TAfterSil1>0)};

    res.dt3sil0=eloss0;
    res.dt3sil1=eloss1;
    res.sildist0=dist0;
    res.sildist1=dist1;
    res.isInSil0=isInSil0;
    res.isInSil1=isInSil1;
    res.silPoint1=silPoint0;

    res.stoppedBeforeSil0=-1;
    res.rangeBeforeSil0=-1;
   /*
    if (isInSil1){
    std::cout<< "silDist1" << silDist1 << "\n";
    }
    */
    return res;
}
