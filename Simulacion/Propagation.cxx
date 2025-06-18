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

    int silIndex0 {-1};
    ROOT::Math::XYZPoint silPoint0 {};
    std::string layer0 {};
    for(const auto& layer : {"f0","l0","r0"})
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

    // Estudiamos el silicio 1: 
    //auto [silIndex0, silPoint0]{sils->FindSPInLayer("f0", vertex, direction)};

    auto silDist0{(vertex - silPoint0).R()};

    double t3inSil0{srim.Slow("light", t3, silDist0)};

    if (incIdx == 0 || incIdx==1){
        t3inSil0=ApplyStraggling(&srim, "light", t3, silDist0);
    }
    

    // Aplicamos straggling y perdidadas de energía a la partícula 3
    auto silThick0{sils->GetLayer(layer0).GetUnit().GetThickness()};
    //std::cout<<"error slow 1 \t"<< silThick0 << " \t t3:" << t3inSil0 <<" \t theta3:" << theta3  << " \n ";

    auto normal {sils->GetLayer(layer0).GetNormal()};
    double thetaSil{TMath::ACos(direction.Unit().Dot(normal.Unit()))};

    auto t3outSil0 {srim.Slow("lightInSil", t3inSil0, silThick0+0.001, thetaSil)};
    double dT3Sil0 {t3inSil0 - t3outSil0};

    if (incIdx == 0 || incIdx==1){
        t3outSil0 =srim.SlowWithStraggling("lightInSil", t3inSil0, silThick0, thetaSil);
        dT3Sil0 = ApplySilResolution(t3inSil0 - t3outSil0);
    }

    // Vemos si la partícula se ha frenaado en el silicio y  si no se ha parado antes
    bool isInSil0{silIndex0 != -1 && 0 < t3inSil0 && t3outSil0 < 0.000001};

    // Vemos si la partícula se ha frenado antes del primer silicio

    double rangeBeforeSil0{(srim.EvalRange("light", t3))};
    bool stoppedBeforeSil0{silDist0 > rangeBeforeSil0};

      
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Estudiamos el silicio 2:
    auto [silIndex1, silPoint1]{sils->FindSPInLayer("f1", vertex, direction)};

    // Calculemos las distancias 
    auto silDist1{(silPoint0 - silPoint1).R()};

    // Aplicamos straggling al t3

    double t3inSil1{srim.Slow("light", t3outSil0, silDist1+0.001)};

    if (incIdx == 0 || incIdx==1){
        t3inSil1=ApplyStraggling(&srim, "light", t3outSil0, silDist1);
    }

    auto silThick1{sils->GetLayer("f1").GetUnit().GetThickness()};
    
    // Aplicamos straggling y perdidadas de energía a la partícula 3


    double t3outSil1{srim.Slow("lightInSil", t3inSil1, silThick1, thetaSil * TMath::DegToRad())};

    // Aplicamos la resolución
    double dT3Sil1{t3inSil1 - t3outSil1};
    if (incIdx == 0 || incIdx==1){
        t3outSil1 =srim.SlowWithStraggling("lightInSil", t3inSil1, silThick1, thetaSil);
        dT3Sil1=ApplySilResolution(t3inSil1 - t3outSil1);
    }


    bool isInSil1{silIndex0 != -1 && 0 < t3inSil0 && 0 < t3inSil1 && silIndex1!=-1 && t3outSil1<0.000001};
    //%%%%%%%%%%%layer0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Guardamos los datos: 

    res.dt3sil0=dT3Sil0;
    res.dt3sil1=dT3Sil1;
    res.sildist0=silDist0;
    res.sildist1=silDist1;
    res.isInSil0=isInSil0;
    res.isInSil1=isInSil1;
    res.silPoint1=silPoint0;

    res.stoppedBeforeSil0=stoppedBeforeSil0;
    res.rangeBeforeSil0=rangeBeforeSil0;
   /*
    if (isInSil1){
    std::cout<< "silDist1" << silDist1 << "\n";
    }
    */
    return res;
}
