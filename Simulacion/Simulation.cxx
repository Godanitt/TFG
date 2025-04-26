#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"


#include "Sampleo.cxx"


//#######################################################
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
 Aplicamos resolucion relacionada con las pérdidas de enerǵia en el gas. Entrada:
 - SRIM::*srim (puntero) -> Nos dice el gas en el que está y como van actuar las pérdidas
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
                                  verticeCentrado.Z() * 10 + 175};

    return verticeUsual;
}

//####################################################### 
/*
Esta es la función principal que se encarga de decir si efectivamente se ha producido la colisión y de obtener histogramas y gráficas
Como argumentos:
 - tBeam: energía del haz (en MeV) (por nucleon)
*/
void Simulation(double tBeam=7.5,double Ex=0.0)
{

    //SetEstiloPublicacionHisto();
    // Definimos las particulas
    Particula li11;
    Particula li10;
    Particula li10ex1;
    Particula li10ex2;
    Particula d;
    Particula t;


    // Detallamos las partículas
    li11.definir_AZ(11, 3,40.79586 ,0.0);
    li10.definir_AZ(10, 3,33.0522, Ex);
    d.definir_AZ(2, 1, 13.13572, 0.0);
    t.definir_AZ(3, 1, 14.94981, 0.0);
    Particula p1{li11};
    Particula p2{d};
    Particula p3{t};
    Particula p4{li10};

    /*
    
    // Definimos la colisión y la detallamos
    Colision colision{li11, d, t, li10,p1.get_A()*tBeam}; 


    // Definimos el detector
    ActarTPC detector;
    double xACTAR{detector.get_xActar()};
    double yACTAR{detector.get_yActar()};
    double zACTAR{detector.get_zActar()};
    std::cout << "Dimensiones del detector: " << xACTAR << " " << yACTAR << " " << zACTAR << '\n';


    // Definimos cinemática de la simulacion 
    auto* Li10aux {new ActPhysics::Kinematics("11Li", "d", "t", p1.get_A() * tBeam)};
    auto* Li10ex1{new ActPhysics::Kinematics("11Li", "d", "t", p1.get_A() * tBeam,5.0)};
    auto* Li10ex2 {new ActPhysics::Kinematics("11Li", "d", "t", p1.get_A() * tBeam,10.0)};
    
    */

    // Grafica de t3 vs theta3 y T4 vs theta4; Sampleamos cinemática, diferentes estados excitados para comprobar

    /*
    SampleoCinematica(li11, d, t, li10, Li10aux, tBeam,0.0);  
    li10ex1.definir_AZ(10, 3,33.0522, 5);
    SampleoCinematica(li11, d, t, li10ex1, Li10ex1, tBeam,5);  
    li10ex2.definir_AZ(10, 3,33.0522, 10);
    SampleoCinematica(li11, d, t, li10ex2, Li10ex2, tBeam,10);  
    */

    /*
     Ahora vamos a reproducir lo que viene siendo la cinemática construida.
    
    Para esto vamos a crear un punto vértica centrado en y,z=0 (gaussiana al rededor de (0,0)) y con x\in(0,xActar).

    Luego vamos a obtener:
     - Espectro de energía (cuentas por energía)
     - Punto de impacto en el silicio 
     - Grafica T vs theta (recostruida)

     Definimos variables de interes antes: 
     - interacciones: número de colisiones 
     - theta3cm: ángulo de salida de la partícula ligera
     - phi3cm: ángulo de salida de la partícula ligera
     - t3: energía de la partícula ligera (salida)
    */
    //#############################################################
    // ACTPHYSICS
    ActRoot::TPCParameters tpc {"Actar"};
    // Tamanhos en mm:
    auto xActar {tpc.X()}; // o mesmo para Y,Z
    auto yActar {tpc.X()}; // o mesmo para Y,Z
    auto zActar {tpc.X()}; // o mesmo para Y,Z
    // Silicios
    auto* sils {new ActPhysics::SilSpecs};
    // Están especificados nun ficheiro de configuración
    sils->ReadFile("./Inputs/silicons_new.conf");
    // Podes comprobar que funciona
    //sils->DrawGeo();

    // Cinemática
    auto* kin {new ActPhysics::Kinematics("11Li", "d", "t", p1.get_A() * tBeam, Ex)}; // cinemática

    // leemos la geometría 
    //auto *geo{new ActSim::Geometry};
    //geo->ReadGeometry("Geometria/", name);


    // SRIM (perdidas de enerǵia)
    ActPhysics::SRIM srim; // &srim -> asi se ''hace'' un pointer
    srim.ReadTable("beam", "SRIM/11Li_mixture_900mbar.txt");
    srim.ReadTable("heavy", "SRIM/10Li_mixture_900mbar.txt");
    srim.ReadTable("heavySil", "SRIM/10Li_silicon.txt");
    srim.ReadTable("light", "SRIM/3H_mixture_900mbar.txt");
    srim.ReadTable("lightInSil", "SRIM/3H_silicon.txt");

    //#############################################################
   

    int interacciones = 1000000;  // numero000 de interacciones
    int fwhm{1}; // FWHM de la resolución
    ROOT::Math::XYZVector vecNormalSil{1, 0, 0}; // vector normal a los siliciones,

    // Histogramas y Gráficas:

    // TH2D: "Nombres","numero bins eje x", "minimo x", "maximo x", "numero bins y", "minimo y", "maximo y"

    auto *histoTheta3T3measured{new TH2D{"histoTheta3T3measured", "T_{4} vs #theta Pre Sil;#theta_{4};T_{4};Contas", 220, 0, 80, 220, -5, 80}};
    auto *histoImpactosSilicioLayer0{new TH2D{"histoImpactosSilicioLayer0", "Impactos Silicio Layer 0; y[mm];z[mm];Contas", 200, -25, 275, 200, 50, 175*2-50}};
    auto *histoImpactosSilicioLayer1{new TH2D{"histoImpactosSilicioLayer1", "Impactos Silicio Layer 1; y[mm];z[mm];Contas", 200, -25, 275, 200, 50, 175*2-50}};
    auto *histoEnergiaParticula3{new TH1D{"histoEnergiaParticula3", "Impactos Silicio; T_{3} [MeV];Cuentas", 200, 0, 80}};
    auto *histoThetaCM3{new TH1D{"histoThetaCM3", "Impactos Silicio; #theta_{CM}[^{#circ}];Cuentas", 200, 0, 70}};
    auto *histoVertex{new TH2D{"histoVertex", "Vertice; z;y", 200, 0, 250,200,0,250}};
    auto *histoDT{new TH2D{"histoDT", "#Delta T_{3} vs T_{3}; T_{3};#Delta T_{3}", 200, 0, 80,200,0,80}};
    auto *histoKinSampled{new TH2D{"histoKinSampled", "#theta_{3lab} vs T_{3lab}; T_{3lab};#theta_{3lab}", 220, 0, 80,220,-5,80}};

    // Saving con Trees 
    auto* outfile {new TFile {TString::Format("./Outputs/tree_Ex_%.2f.root", Ex), "recreate"}};
    auto* tree {new TTree {"SimulationTree", "Simple simulation tree"}};
    double t3rec_tree {};
    tree->Branch("t3rec", &t3rec_tree);
    double exrec_tree {};
    tree->Branch("exrec", &exrec_tree);
    double theta3_tree {};
    tree->Branch("theta3", &theta3_tree);

    for (int i=0; i<interacciones; i++){

        // definimos el vértice 
        ROOT::Math::XYZPoint vertex{SampleVertex(xActar, yActar, zActar)};
        //std::cout << "Vetice:" << vertex.X() << " " << vertex.Y() << " " << vertex.Z() << '\n';
        histoVertex->Fill(vertex.X(), vertex.Y());
        // Reducimos la energía. 1-> Calulamos la posición del impacto. 2-> Nueva enerǵia del haz.
        auto r{vertex.X()}; 
        auto tBeamNew{srim.Slow("beam", tBeam*p1.get_A(), r, 0 * TMath::DegToRad())};
        //std::cout << "TbeamNew: " << tBeamNew << '\n';


        // Calculamos el ángulo theta3 y phi3 cm. Para esto creamos los ángulos aleatoriamente 
        auto uniforme1{gRandom->Uniform(-1,1)};
        auto uniforme2{gRandom->Uniform(0,2*TMath::Pi())};
        // Aquí iría la seccion eficaz, de momento d#sigma vs d#Omega = 1 
        auto theta1CM{TMath::ACos(uniforme1)};
        auto phiCM{uniforme2};
        auto [t3,theta3] {valores3(p1, p2, p3, p4, tBeamNew,theta1CM )};

        histoKinSampled->Fill(theta3*TMath::RadToDeg(), t3);

        theta3 = ApplyAngleResolution(theta3, fwhm*TMath::DegToRad());

        

        //std::cout << "Theta3: " << theta3* TMath::RadToDeg()<< " T3 : " << t3 << '\n';


        // Proyectamos la posición del ángulo incidente a los detectores: 
        ROOT::Math::XYZVector direction{std::cos(theta3), std::sin(theta3)*std::sin(phiCM), std::sin(theta3)*std::cos(phiCM)};
      
        // Para comprobar que a partícula chega a un silicio:
        // 1o argumento: nome da layer á que queres ver se chegas. a primeira layer é "f0"
        // 2o vértice
        // 3o dirección
        // Para a dirección ten en conta que en ACTAR TPC traballamos con X <--> Z.
        // O eixo de propagación do beam é X
        // Lembra que a dirección é o vector unitario r en coordenadas ESFÉRICAS
        // Se silIndex == -1 -> NON HOUBO IMPACTO
        // Se silIndex != 1 -> SI houbo impacto, no punto silPoint0
        auto [silIndex0, silPoint0] {sils->FindSPInLayer("f0", vertex, direction)};

        //Calculamos el grosor del silicio
        auto silDist0 {(vertex - silPoint0).R()};
        /*
        double dx = vertex.X() - silPoint0.X(); 
        double dy = vertex.Y() - silPoint0.Y();
        double dz = vertex.Z() - silPoint0.Z();
        double silDist = std::sqrt(dx*dx + dy*dy + dz*dz);
        */
        if (silIndex0 != -1)
        {
            // Guardamos las gráficas de las posiciones en

            histoThetaCM3->Fill(theta3 * TMath::RadToDeg());
            histoImpactosSilicioLayer0->Fill(silPoint0.Y(), silPoint0.Z());

            // Calculamos el ángulo de incidencia de 4 respecto al silicio
            double thetaSi{TMath::ACos(direction.Dot(vecNormalSil))};

            // Aplicamos straggling al t3
            double t3inSil0{ApplyStraggling(&srim, "light", t3, silDist0)};

            // Check particle is not stopped in gas
            if (t3inSil0 < 0)
                continue;

            // Aplicamos straggling y perdidadas de energía a la partícula 3
            auto silThick0 {sils->GetLayer("f0").GetUnit().GetThickness()};  
            double t3outSil0{srim.Slow("lightInSil", t3inSil0, silThick0, thetaSi* TMath::DegToRad() )};
            double dT3Sil0{ApplySilResolution(t3inSil0 - t3outSil0)};

            // True -> Es frenada por el primer detector (silicio 0), False -> NO es frenada  -> Recostruimos la cinemática incluyendo silicio 1:
            bool isInSil{t3outSil0<0.0001};
            std::cout << "t3outSil0: " << t3outSil0 << '\n';

            // Calculamos la energían inicial en el vértice según nuestro modelo
            double t3Measured0{srim.EvalInitialEnergy("light", dT3Sil0, silDist0)};

            // Si es frenada -> Recostruimos cinemática 
            if (isInSil)
            {
                std::cout << "cayo en el primero ;(): " << '\n';
                // Guardamos energías
                histoEnergiaParticula3->Fill(t3Measured0);
                histoTheta3T3measured->Fill(theta3 * TMath::RadToDeg(), t3Measured0);
                histoDT->Fill(t3Measured0, dT3Sil0);
                
                //Reconstruction
                auto recEx {kin->ReconstructExcitationEnergy(t3Measured0, theta3)};
                auto recThetaCM {kin->ReconstructTheta3CMFromLab(t3Measured0, theta3)};
    
                // Guardamos
                t3rec_tree = t3Measured0;
                exrec_tree = recEx;
                theta3_tree = theta3;
                tree->Fill();
            }
            else {
                /* 
                Segundo detector. Vamos a evaluar si llega al segundo detector para recostruir correctamente la energía de la partícula saliente:         
                */
                auto [silIndex1, silPoint1] {sils->FindSPInLayer("f1", vertex, direction)};
                if (silIndex1 != -1)
                {
                    histoImpactosSilicioLayer1->Fill(silPoint1.Y(), silPoint1.Z());
                    // Calculamos la distancia entre el punto en el que colisiona en el primer silicio y en el segundo
                    auto silDist1 {(silPoint0 - silPoint1).R()};

                    // Aplicamos straggling al t3
                    double t3inSil1{ApplyStraggling(&srim, "light", t3outSil0, silDist1)};
        
                    // Comprobamos si la partícula no se para en el InterSil
                    if (t3inSil1 < 0){
                        continue;
                    }
                    // Aplicamos straggling y perdidadas de energía a la partícula 3
                    auto silThick1 {sils->GetLayer("f1").GetUnit().GetThickness()};  
                    double t3outSil1=srim.Slow("lightInSil", t3inSil1, silThick1, thetaSi* TMath::DegToRad());
                    // Aplicamos la resolución
                    double dT3Sil1{ApplySilResolution(t3inSil1 - t3outSil1)};
        
                    //Calculamos el grosor del silicio
                    // Calculamos la energían inicial en el vértice según nuestro modelo
                    double t3Measured1{srim.EvalInitialEnergy("light", dT3Sil1, silDist1)};
                    std ::cout << "t3Measured0: " << t3Measured0 << '\n';
                    std ::cout << "t3Measured1: " << t3Measured1 << '\n';
                    t3Measured1+=t3Measured0;
                    std ::cout << "t3Measured2: " << t3Measured1 << '\n';
    
                    // Guardamos energías
                    histoEnergiaParticula3->Fill(t3Measured1);  
                    histoTheta3T3measured->Fill(theta3 * TMath::RadToDeg(), t3Measured1);
                    histoDT->Fill(t3Measured1, dT3Sil0);
 
                    //Reconstruction
                    auto recEx {kin->ReconstructExcitationEnergy(t3Measured1, theta3)};
                    auto recThetaCM {kin->ReconstructTheta3CMFromLab(t3Measured1, theta3)};
        

                    // Guardamos
                    t3rec_tree = t3Measured1;
                    exrec_tree = recEx;
                    theta3_tree = theta3;
                    tree->Fill();

                }
            }
        }
}

outfile->cd();
tree->Write();

// Graficamos (temporal)
gStyle->SetPadRightMargin(0.15);  // margen derecho
// Graficamos el T3 vs theta3 measured
auto *c1{new TCanvas{"c1", "T3_vs_theta"}};


c1->DivideSquare(6);
c1->cd(1);
auto* g3Li10 {kin->GetKinematicLine3()};
histoTheta3T3measured->Draw("colz");
g3Li10->Draw("l");
c1->cd(1)->SaveAs("Graficas/Cinematica Reconstruida/CinematicaSimulation.eps");
// Vertice: 
c1->cd(2);
histoVertex->Draw("colz");
c1->cd(2)->SaveAs("Graficas/Cinematica Reconstruida/VerticeXY.eps");
// Graficamos punto impacto en el silicio ;
c1->cd(3);
histoImpactosSilicioLayer0->Draw("colz");
sils->GetLayer("f0").GetSilMatrix()->Draw();
c1->cd(3)->SaveAs("Graficas/Impacto/Impacto_Layers1.eps");
// Graficamos punto impacto en el silicio ;
c1->cd(4);
histoImpactosSilicioLayer1->Draw("colz");
sils->GetLayer("f1").GetSilMatrix()->Draw();
c1->cd(4)->SaveAs("Graficas/Impacto/Impacto_Layers2.eps");
// Graficamos la perdida de energía
c1->cd(5);
histoDT->Draw("colz");
c1->cd(5)->SaveAs("Graficas/Cinematica Reconstruida/Perdidas_Energia.eps");
// Graficamos la cinemática sampleada;
c1->cd(6);
histoKinSampled->Draw("colz");
c1->cd(6)->SaveAs("Graficas/Cinematica Sampleada/CinematicaSampleada.eps");



}   
//####################################################### 