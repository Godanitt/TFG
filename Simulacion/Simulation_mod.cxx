#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"

#include "TEfficiency.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

#include "Auxiliares.cxx"
#include "Propagation.cxx"
#include "Sampleo.cxx"

/*
Esta es la función principal que se encarga de decir si efectivamente se ha 
roducido la colisión y de obtener histogramas y gráficas Como argumentos:
 - tBeam: energía del haz (en MeV) (por nucleon)
*/
void Simulation_mod(double tBeam = 7.5, double Ex = 0.0, int interacciones = 1000000,int incIdx=0)
{
    // Definimos las particulas
    Particula li11;
    Particula li10;
    Particula li10ex1;
    Particula li10ex2;
    Particula d;
    Particula t;

    // Detallamos las partículas
    li11.definir_AZ(11, 3, 40.79586, 0.0);
    li10.definir_AZ(10, 3, 33.0522, Ex);
    d.definir_AZ(2, 1, 13.13572, 0.0);
    t.definir_AZ(3, 1, 14.94981, 0.0);
    Particula p1 {li11};
    Particula p2 {d};
    Particula p3 {t};
    Particula p4 {li10};

    /*
     Ahora vamos a reproducir lo que viene siendo la cinemática construida.

    Para esto vamos a crear un punto vértica centrado en y,z=0 (gaussiana al 
    ededor de (0,0)) y con x\in(0,xActar).

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

    // #############################################################

    // ACTPHYSICS

    // Actar
    ActRoot::TPCParameters tpc {"Actar"};
    // Tamanhos en mm:
    auto xActar {tpc.X()}; // o mesmo para Y,Z
    auto yActar {tpc.Y()}; // o mesmo para Y,Z
    auto zActar {tpc.Z()}; // o mesmo para Y,Z

    // Silicios
    auto* sils {new ActPhysics::SilSpecs};
    sils->ReadFile("./Inputs/silicons_new.conf");
    sils->GetLayer("f0").MoveZTo(75, {3});
    sils->GetLayer("f1").MoveZTo(75, {3});
    sils->GetLayer("l0").MoveZTo(75, {3});
    sils->GetLayer("r0").MoveZTo(75, {3});

    sils->DrawGeo();

    // Cinemática
    auto* kin {new ActPhysics::Kinematics("11Li", "d", "t", p1.get_A() * tBeam, Ex)}; // cinemática

    // SRIM (perdidas de enerǵia)
    ActPhysics::SRIM srim; // &srim -> asi se ''hace'' un pointer
    srim.SetUseSpline(false);
    srim.ReadTable("beam", "SRIM/11Li_mixture_900mbar.txt");
    srim.ReadTable("heavy", "SRIM/10Li_mixture_900mbar.txt");
    srim.ReadTable("heavySil", "SRIM/10Li_silicon.txt");
    srim.ReadTable("light", "SRIM/3H_mixture_900mbar.txt");
    srim.ReadTable("lightInSil", "SRIM/3H_silicon.txt");

    // #############################################################

    int fwhm {1};                                 // FWHM de la resolución
    ROOT::Math::XYZVector vecNormalSil {1, 0, 0}; // vector normal a los siliciones,

    // Histogramas y Gráficas:
    // TH2D: "Nombres","numero bins eje x", "minimo x", "maximo x", "numero bins
    // y", "minimo y", "maximo y"
    auto* hRPz {new TH2D {"hRPz", "RP;X;Z", 350, 0, 270, 350, 0, 270}};
    auto* hKinSampled {
        new TH2D {"hKinSampled", "#theta_{3} vs T_{3} Sampled;#theta_{3}; T_{3}", 220, 0, 80, 220, -5, 80}};
    auto* hKinMeasured {
        new TH2D {"hKinMeasured", "#theta_{3} vs T_{3} Recostrued; #theta_{3};T_{3};Contas", 220, 0, 80, 220, -5, 80}};
    auto* hImpactSilF0 {new TH2D {"hImpactSilF0", "Impactos Silicio Layer f0; y[mm];z[mm];Contas", 200, -25, 275, 200,
                                  0, 175 * 2 - 100}};
    auto* hImpactSilF1 {new TH2D {"hImpactSilF1", "Impactos Silicio Layer f1; y[mm];z[mm];Contas", 200, -25, 275, 200,
                                  0, 175 * 2 - 100}};
    auto* hImpactSilL0 {new TH2D {"hImpactSilL0", "Impactos Silicio Layer l0; y[mm];z[mm];Contas", 250, -200, 400, 200,
                                  -50, 175 * 2 - 100}};
    auto* hImpactSilR0 {new TH2D {"hImpactSilR0", "Impactos Silicio Layer r0; y[mm];z[mm];Contas", 250, -200, 400, 200,
                                  -50, 175 * 2 - 100}};
    // TH1D
    auto* hThetaCMAll {new TH1D {"hThetaCMAll", "All #theta_{CM}; #theta; Contas", 180, 0, 180}};
    auto* hThetaCMIn {new TH1D {"hThetaCMIn", "#theta_{CM} in cuts; #theta; Contas", 180, 0, 180}};
    auto* hThetaCMRec {new TH1D {"hThetaCMRec", "#theta_{CM} rec;  #theta; Contas", 180, 0, 180}};
    auto* hEx {new TH1D {"hEx", "Rec Ex;E_{x} [MeV];Counts", 300, -5, 5}};

    // Saving con Trees
    auto* outfile {new TFile {TString::Format("./Outputs/tree_Ex_%.2f.root", Ex), "recreate"}};
    auto* tree {new TTree {"SimulationTree", "Simple simulation tree"}};
    double t3rec_tree {};
    tree->Branch("t3rec", &t3rec_tree);
    double exrec_tree {};
    tree->Branch("exrec", &exrec_tree);
    double theta3_tree {};
    tree->Branch("theta3", &theta3_tree);

    // Simulacion core
    for(int i = 0; i < interacciones; i++)
    {

        // Definimos el vértice
        ROOT::Math::XYZPoint vertex {SampleVertex(xActar, yActar, zActar)};
        hRPz->Fill(vertex.X(), vertex.Z());

        // Reducimos la energía. 1-> Calulamos la posición del impacto. 2-> Nueva e
        // erǵia del haz.
        auto r {vertex.X()};
        auto tBeamNew {srim.Slow("beam", tBeam * p1.get_A(), r, 0 * TMath::DegToRad())};

        // Calculamos el ángulo theta3 y phi3 cm. Para esto creamos los ángulos a
        // eatoriamente
        auto uniforme1 {gRandom->Uniform(-1, 1)};
        auto uniforme2 {gRandom->Uniform(0, 2 * TMath::Pi())};

        // Aquí iría la seccion eficaz, de momento d#sigma vs d#Omega = 1 
        auto thetaCM {TMath::ACos(uniforme1)};
        hThetaCMAll->Fill(thetaCM * TMath::RadToDeg());
        auto phiCM {uniforme2};

        // Calculamos la cinemática
        kin->ComputeRecoilKinematics(thetaCM, phiCM);
        double t3 {kin->GetT3Lab()};
        double theta3 {kin->GetTheta3Lab()};

        // Aplicamos resolución angular
        //theta3 = ApplyAngleResolution(theta3, fwhm * TMath::DegToRad());

        // Proyectamos la posición del ángulo incidente a los detectores: 200
        ROOT::Math::XYZVector direction {std::cos(theta3), std::sin(theta3) * std::sin(phiCM),
                                         std::sin(theta3) * std::cos(phiCM)};
        hKinSampled->Fill(theta3 * TMath::RadToDeg(), t3);

        // Obtenemos los silIndex
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
        if(silIndex0 == -1 || layer0 == "")
            continue;
        // Propagate to it
        auto dist0 {(vertex - silPoint0).R()};
        auto TAtSil0 {srim.Slow("light", t3, dist0)};
        // If stops before, skip event too
        if(TAtSil0 <= 0)
            continue;
        ///////////// SILICON 0 ///////////////////
        // Angle with normal
        auto thick0 {sils->GetLayer(layer0).GetUnit().GetThickness()};
        auto normal {sils->GetLayer(layer0).GetNormal()};
        auto angleNormal0 {TMath::ACos(direction.Unit().Dot(normal.Unit()))};
        // Eloss
        auto TAfterSil0 {srim.Slow("lightInSil", TAtSil0, thick0, angleNormal0)};
        auto eloss0 {TAtSil0 - TAfterSil0};
        // std::cout << "TAtSil0 : " << TAtSil0 << '\n';
        // std::cout << "thick0 : " << thick0 << '\n';
        // std::cout << "angleNormal0 : " << angleNormal0 << '\n';
        // std::cout << "eloss0 : " << eloss0 << '\n';

        //////////// SILICON 1 ////////////////////
        // only if layer0 == "f0"
        bool isInSil1 {};
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
            if(silIndex1 == -1)
                continue;
            // Eloss in intergas = region between f0 and f1
            dist1 = (silPoint0 - silPoint1).R();
            TAtSil1 = srim.Slow("light", TAfterSil0, dist1);
            if(TAtSil1 <= 0) // event undergoes punch in f0 but stops before f1
                continue;
            auto thick1 {sils->GetLayer("f1").GetUnit().GetThickness()};
            auto normal1 {sils->GetLayer("f1").GetNormal()};
            auto angleNormal1 {TMath::ACos(direction.Unit().Dot(normal.Unit()))};
            TAfterSil1 = srim.Slow("lightInSil", TAtSil1, thick1, angleNormal1);
            eloss1 = TAtSil1 - TAfterSil1;
            // Punch outside sil1
            if(TAfterSil1 > 0)
                continue;
            isInSil1 = true;
        }
        if(layer0 != "f0" && TAfterSil0 > 0) // skip punch for layers other than f0
            continue;

        // Reconstruct
        double recTAtSil0 {};
        if(isInSil1)
        {
            auto recTAfterSil0 {srim.EvalInitialEnergy("light", eloss1, dist1)};
            recTAtSil0 = eloss0 + recTAfterSil0;
        }
        else
            recTAtSil0 = eloss0;
        auto recT3 {srim.EvalInitialEnergy("light", recTAtSil0, dist0)};
        auto recEx {kin->ReconstructExcitationEnergy(recT3, theta3)};
        auto recThetaCM {kin->ReconstructTheta3CMFromLab(recT3, theta3)};
        // Fill histograms
        hKinMeasured->Fill(theta3 * TMath::RadToDeg(), recT3);
        hThetaCMIn->Fill(thetaCM * TMath::RadToDeg());
        hEx->Fill(recEx);

        if(layer0 == "f0")
            hImpactSilF0->Fill(silPoint0.Y(), silPoint0.Z());
        if(layer0 == "f0" && isInSil1)
            hImpactSilF1->Fill(silPoint0.Y(), silPoint0.Z());
        if(layer0 == "l0")
            hImpactSilL0->Fill(silPoint0.X(), silPoint0.Z());
        if(layer0 == "r0")
            hImpactSilR0->Fill(silPoint0.X(), silPoint0.Z());

        // auto [silIndex1, silPoint1] {sils->FindSPInLayer("f1", vertex, direction)};
        //
        // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // // Propagamos (véase Propagation.cxx)
        // Resultado propagacion = Propagation(t3, theta3, vertex, direction, srim, sils);
        //
        // // Perdidas de energia
        // double dT3Sil0 = propagacion.dt3sil0;
        // double dT3Sil1 = propagacion.dt3sil1;
        //
        // // Distancias en las que se perdieron la energía:
        // double silDist0 = propagacion.sildist0;
        // double silDist1 = propagacion.sildist1;
        //
        // // Puntos de imapacto del silicio:
        // // ROOT::Math::XYZPoint silPoint1=propagacion.silPoint1;
        //
        // // Booleano que nos dice si se ha medido en el layer 0:
        // bool isInSil0 = propagacion.isInSil0;
        //
        // // Booleano que nos dice si se ha medido en el layer 1:
        // bool isInSil1 = propagacion.isInSil1;
        //
        // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // // Recostruimos:

        // bool isInSil0 {TAfterSil0 <= 0};
        // bool isInSil1 {TAtSil1 > 0 && TAfterSil1 <= 0};
        //
        // double recT3;
        // double recEx;
        // double recThetaCM;
        //
        // if(isInSil0)
        // {
        //     recT3 = srim.EvalInitialEnergy("light", eloss0, dist0);
        //     recEx = kin->ReconstructExcitationEnergy(recT3, theta3);
        //     recThetaCM = kin->ReconstructTheta3CMFromLab(recT3, theta3);
        // }
        //
        // if(isInSil1)
        // {
        //     std::cout << "insil1" << '\n';
        //     // Order matters
        //     // 1-> Particle has stopped on sil1. Therefore, its total energy is dT3Sil1
        //     // We need energy after sil0
        //     auto recTAfterSil0 {srim.EvalInitialEnergy("light", eloss1, dist1)};
        //     // 2-> But then energy to recover is eloss0 + recTAfterSil0!
        //     auto recT3AtSil0 {dist0 + recTAfterSil0};
        //     // 3-> And this is the total energy to recover at vertex
        //     // double recT3sil0 {srim.EvalInitialEnergy("light", dT3Sil0, silDist0)};
        //     // double recT3sil1 {srim.EvalInitialEnergy("light", dT3Sil1, silDist1)};
        //     // recT3 = recT3sil0 + recT3sil1;
        //     recT3 = srim.EvalInitialEnergy("light", recT3AtSil0, dist0);
        //
        // }
        //
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Enchemos


        // if(isInSil0 || isInSil1)
        // {
        //     // Rellenamos:
        //     hKinMeasured->Fill(theta3 * TMath::RadToDeg(), recT3);
        //     hThetaCMRec->Fill(recThetaCM);
        //
        //     // Rellenamos tree
        //     t3rec_tree = recT3;
        //     exrec_tree = recEx;
        //     theta3_tree = theta3;
        //     tree->Fill();
        // }
    }

    outfile->cd();
    tree->Write();

    // Graficamos el T3 vs theta3 measured
    auto* c1 {new TCanvas {"c1", "T3_vs_theta"}};
    gStyle->SetPadRightMargin(0.15);

    c1->DivideSquare(4);

    c1->cd(1);
    hKinMeasured->Draw("colz");

    c1->cd(2);
    hKinSampled->Draw("colz");

    c1->cd(3);
    auto eff {new TEfficiency(*hThetaCMIn, *hThetaCMAll)};
    eff->SetTitle("Eficiencia;#theta_{CM} [#circ];#epsilon");
    eff->Draw("AP");
    c1->cd(4);
    hEx->Draw();

    // c1->SaveAs(TString::Format("Graficas/Completo_Ex%.2f.pdf", Ex));
    // c1->SaveAs(TString::Format("Graficas/Completo_Ex%.2f.eps", Ex));

    // eff->GetPaintedGraph()->GetXaxis()->SetTitle("#theta_{CM}");
    // eff->GetPaintedGraph()->GetYaxis()->SetTitle("Eficiencia");

    // Dibujamos Angulos

    /*

    auto *c2{new TCanvas{"c2", "theta cm"}};
    c2->cd();
    static int colorIndex = 2;  // Empieza en 2 (rojo)
    hThetaCMSampled->Draw();
    hThetaCMSampled->SetFillColor(colorIndex);
    hThetaCMRec->Draw("same");
    hThetaCMRec->SetFillColor(colorIndex+1);
    c2->Update();
    c2->SaveAs(TString::Format("Graficas/ThetaCM_Ex%.2f.pdf", Ex));
    c2->SaveAs(TString::Format("Graficas/ThetaCM_Ex%.2f.eps", Ex));
    */

    // Dibujamos puntos de impacto
    auto* c3 {new TCanvas {"c3", "T3_vs_theta"}};
    c3->DivideSquare(4);
    c3->cd(1);
    hImpactSilF0->Draw("colz");
    sils->GetLayer("f0").GetSilMatrix()->Draw();

    c3->cd(2);
    hImpactSilF1->Draw("colz");
    sils->GetLayer("f1").GetSilMatrix()->Draw();

    c3->cd(3);
    hImpactSilL0->Draw("colz");
    sils->GetLayer("l0").GetSilMatrix()->Draw();

    c3->cd(4);
    hImpactSilR0->Draw("colz");
    sils->GetLayer("r0").GetSilMatrix()->Draw();

    // c3->SaveAs(TString::Format("Graficas/Impacts_Ex%.2f.pdf", Ex));
    // c3->SaveAs(TString::Format("Graficas/Impacts_Ex%.2f.eps", Ex));
}
// #######################################################
