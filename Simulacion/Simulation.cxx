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

#include "TVector3.h" 
#include "Math/Point3D.h"
#include "Math/Vector3D.h"

#include "Sampleo.cxx"
#include "Propagation.cxx"
#include "Auxiliares.cxx"

/*
Esta es la función principal que se encarga de decir si efectivamente se ha producido la colisión y de obtener histogramas y gráficas
Como argumentos:
 - tBeam: energía del haz (en MeV) (por nucleon)
 - Ex: energía de excitacion
 - interacciones: numero de particulas en el haz 
 - indc =  Array de índices para seleccionar diferentes tipos de incertidumbre. Cada elemento representa un identificador del tipo de incertidumbre:
    0: todas las fuentes de incertidumbre
    1: Solo straggling
    2: Solo theta
    3: ninguna fuente de incertidumbre
 */

void Simulation(double tBeam=7.5,double Ex=0.0, int interacciones = 1000000, int incIdx = 0)
{
    double ExOg {Ex};
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

    // Actar
    ActRoot::TPCParameters tpc {"Actar"};
    // Tamanhos en mm:
    auto xActar {tpc.X()}; // o mesmo para Y,Z
    auto yActar {tpc.X()}; // o mesmo para Y,Z
    auto zActar {tpc.X()}; // o mesmo para Y,Z

    // Silicios
    auto* sils {new ActPhysics::SilSpecs};
    sils->ReadFile("./Inputs/silicons_new.conf");
    sils->GetLayer("f0").MoveZTo(75, {3});
    sils->GetLayer("f1").MoveZTo(75, {3});
    sils->GetLayer("l0").MoveZTo(75, {3});
    sils->GetLayer("r0").MoveZTo(75, {3});

    sils->DrawGeo();


    // SRIM (perdidas de enerǵia)
    ActPhysics::SRIM srim; // &srim -> asi se ''hace'' un pointer
    srim.ReadTable("beam", "SRIM/11Li_900mb_CF4_90-10.txt");
    srim.ReadTable("heavy", "SRIM/10Li_900mb_CF4_90-10.txt");
    srim.ReadTable("heavySil", "SRIM/10Li_silicon.txt");
    srim.ReadTable("light", "SRIM/3H_900mb_CF4_90-10.txt");
    srim.ReadTable("lightInSil", "SRIM/3H_silicon.txt");

    //#############################################################

    int fwhm{1}; // FWHM de la resolución
    ROOT::Math::XYZVector vecNormalSil{1, 0, 0}; // vector normal a los siliciones,

    // Histogramas y Gráficas:
    // TH2D: "Nombres","numero bins eje x", "minimo x", "maximo x", "numero bins y", "minimo y", "maximo y"
    auto *hKinSampled{new TH2D{"hKinSampled", "#theta_{3} vs T_{3} Sampled;#theta_{3}; T_{3}", 220, 0, 80,220,-5,80}};
    auto *hKinMeasured{new TH2D{"hKinMeasured", "#theta_{3} vs T_{3} Recostrued; #theta_{3};T_{3};Contas", 220, 0, 80, 220, -5, 80}};
    auto *hImpactSilF0{new TH2D{"hImpactSilF0", "Impactos Silicio Layer f0; y[mm];z[mm];Contas", 200, -25, 275, 200, 0, 175*2-100}};
    auto *hImpactSilF1{new TH2D{"hImpactSilF1", "Impactos Silicio Layer f1; y[mm];z[mm];Contas", 200, -25, 275, 200, 0, 175*2-100}};
    auto *hImpactSilL0{new TH2D{"hImpactSilL0", "Impactos Silicio Layer l0; y[mm];z[mm];Contas", 250, -200, 400, 200, -50, 175*2-100}};
    auto *hImpactSilR0{new TH2D{"hImpactSilR0", "Impactos Silicio Layer r0; y[mm];z[mm];Contas", 250, -200, 400, 200, -50, 175*2-100}};
    // TH1D
    auto *hThetaCMSampled{new TH1D{"hThetaCMSampled", "#theta_{CM} sampled; #theta; Contas", 180, 0, 180}};
    auto *hThetaCMRec{new TH1D{"hThetaCMRec", "#theta_{CM} rec;  #theta; Contas", 180, 0, 180}};
    auto* hEx {new TH1D {"hEx", "Rec Ex;E_{x} [MeV];Counts", 300, -5, 5}};
    auto* hThetaCM {new TH1D {"hThetaCM", "CM;#theta_{CM};Counts", 300, 0, 180}};
    auto* hRange {new TH1D {"hRangue", "Rangue;#Rangue [mm];Counts", 300, 0, 255}};


    // Saving con Trees 
    auto* outfile {new TFile {TString::Format("./Outputs/tree_Ex_%.2f_incIdx%i.root", ExOg,incIdx), "recreate"}};
    auto* tree {new TTree {"SimulationTree", "Simple simulation tree"}};
    double t3rec_tree {};
    tree->Branch("t3rec", &t3rec_tree);
    double exrec_tree {};
    tree->Branch("exrec", &exrec_tree);
    double theta3_tree {};
    tree->Branch("theta3", &theta3_tree);


    // Cinemática
    auto* kin {new ActPhysics::Kinematics("11Li", "d", "t", p1.get_A() * tBeam, Ex)}; // cinemática

    // Anchura de la Briet-Wigner: 
    double Gamma{0.01}; // en MeV
    ExOg=Ex;

    // CrossSection:
    ActSim::CrossSection xs {};
    
    if (Ex==0.0){
        xs.ReadFile("Inputs/xs/s12_p1i.dat");
    }
    if (Ex==0.2){
        xs.ReadFile("Inputs/xs/p12_p1i.dat");
    }


    // Simulacion core
    for (int i=0; i<interacciones; i++){
        if (ExOg==0.0) {
            Gamma=0.1;
            Ex=gRandom->BreitWigner(ExOg, Gamma);
            if (Ex<ExOg-5 || Ex>ExOg+5){
                continue;
            }
            kin=new ActPhysics::Kinematics("11Li", "d", "t", p1.get_A() * tBeam, Ex); // cinemática
        }

        if (ExOg==0.2) {
            Gamma=0.2;
            Ex=gRandom->BreitWigner(ExOg, Gamma);
            if (Ex<ExOg-5 || Ex>ExOg+5){
                continue;
            }
            kin=new ActPhysics::Kinematics("11Li", "d", "t", p1.get_A() * tBeam, Ex); // cinemática
        }
        

        // Definimos el vértice
        ROOT::Math::XYZPoint vertex{SampleVertex(xActar, yActar, zActar)};

        // Reducimos la energía. 1-> Calulamos la posición del impacto. 2-> Nueva enerǵia del haz.
        auto r{vertex.X()}; 
        auto tBeamNew{srim.Slow("beam", tBeam*p1.get_A(), r, 0 * TMath::DegToRad())};

        // Calculamos el ángulo theta3 y phi3 cm. Para esto creamos los ángulos aleatoriamente 
        auto uniforme1{gRandom->Uniform(-1,1)};
        auto uniforme2{gRandom->Uniform(0,2*TMath::Pi())};

        // Aquí iría la seccion eficaz, de momento d#sigma vs d#Omega = 1 
        //auto thetaCM{TMath::ACos(uniforme1)};

        
        double thetaCM{(xs.SampleHist())*TMath::DegToRad()};


        auto phiCM{uniforme2};


        // Calculamos la cinemática
        kin->ComputeRecoilKinematics(thetaCM, phiCM);
        double t3 {kin->GetT3Lab()};
        double theta3 {kin->GetTheta3Lab()};

        // Aplicamos resolución angular
        if (incIdx == 0 || incIdx == 2) {
            theta3 = ApplyAngleResolution(theta3, fwhm*TMath::DegToRad());      
        }

        // Proyectamos la posición del ángulo incidente a los detectores: 200
        ROOT::Math::XYZVector direction{std::cos(theta3), std::sin(theta3)*std::sin(phiCM), std::sin(theta3)*std::cos(phiCM)};

        // Obtenemos los silIndex
        int silIndex0 {};
        ROOT::Math::XYZPoint silPoint0 {};
        std::string layer0 {};
        std::string l0 {"l0"};
        std::string f0 {"f0"};
        std::string r0 {"r0"};
        for(const auto& layer : {"f0","l0","r0"})
        {
            auto [idx, point] = sils->FindSPInLayer(layer, vertex, direction);
            if(idx != -1)
            {
                silIndex0 = idx;
                silPoint0 = point;
                layer0 = layer;
            }
        }
        
        auto [silIndex1, silPoint1]{sils->FindSPInLayer("f1", vertex, direction)};


        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Propagamos (véase Propagation.cxx)
        
        Resultado propagacion = Propagation(t3, theta3, vertex, direction, srim, sils,incIdx);


        // Perdidas de energia
        double dT3Sil0=propagacion.dt3sil0;
        double dT3Sil1=propagacion.dt3sil1;


        // Distancias en las que se perdieron la energía:
        double silDist0=propagacion.sildist0;
        double silDist1=propagacion.sildist1;
        //std::cout << silDist1 << "\n";

        // Puntos de imapacto del silicio:
        //ROOT::Math::XYZPoint silPoint1=propagacion.silPoint1;

        // Booleano que nos dice si se ha medido en el layer 0:
        bool isInSil0=propagacion.isInSil0;

        // Booleano que nos dice si se ha medido en el layer 1:
        bool isInSil1=propagacion.isInSil1;

        // Valores que nos dicen si se ha parado antes de chocar con el silicio + nos da su rango dentro d eeste
        bool stoppedBeforeSil0=propagacion.stoppedBeforeSil0;
        double rangeBeforeSil0=propagacion.rangeBeforeSil0;

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Recostruimos:

        double recT3;
        double recEx;
        double recThetaCM{0.0};


        if (isInSil0) {
            recT3=srim.EvalInitialEnergy("light", dT3Sil0, silDist0);
            recEx =kin->ReconstructExcitationEnergy(recT3, theta3);
            recThetaCM =kin->ReconstructTheta3CMFromLab(recT3, theta3);       
        }
        
        if (isInSil1) {
            double recT3sil0{srim.EvalInitialEnergy("light", dT3Sil0, silDist0)};
            double recT3sil1{srim.EvalInitialEnergy("light", dT3Sil1, silDist1)};

            recT3=recT3sil0+recT3sil1;

            recEx =kin->ReconstructExcitationEnergy(recT3, theta3);
            recThetaCM =kin->ReconstructTheta3CMFromLab(recT3, theta3);
        }      

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Enchemos 

        hThetaCMSampled->Fill(thetaCM*TMath::RadToDeg());
        hKinSampled->Fill(theta3*TMath::RadToDeg(),t3);

        if (isInSil0 || isInSil1) {
            // Rellenamos: 
            hKinMeasured->Fill(theta3*TMath::RadToDeg(),recT3);
            hThetaCMRec->Fill(thetaCM*TMath::RadToDeg());
            hEx->Fill(recEx);

            // Rellenamos tree
            t3rec_tree = recT3;
            exrec_tree = recEx;
            theta3_tree = theta3;
            tree->Fill();     

        }
        

        if (stoppedBeforeSil0){
            hRange->Fill(rangeBeforeSil0);
            //std::cout<<"Range"<<silDist0<<"\n";
            //std::cout<<"Range"<<rangeBeforeSil0<<"\n";
        }


        if (silIndex0!=-1 && layer0==f0){
            hImpactSilF0->Fill(silPoint0.Y(),silPoint0.Z());
        }
        
        if (silIndex0!=-1 && silIndex1!=-1 && layer0==f0){
            hImpactSilF1->Fill(silPoint1.Y(),silPoint1.Z());
        }
        
        if (silIndex0!=-1 && layer0==l0){
            //std::cout << "Silpoit l0 \t" << silPoint0 << "\n";
            hImpactSilL0->Fill(silPoint0.X(),silPoint0.Z());
        }            
        
        if (silIndex0!=-1 && layer0==r0){
            //std::cout << "Silpoit r0 \t" << silPoint0 << "\n";
            hImpactSilR0->Fill(silPoint0.X(),silPoint0.Z());
        }
    }

    outfile->cd();
    tree->Write();


    // Graficamos el T3 vs theta3 measured
    auto *c1{new TCanvas{"c1", "T3_vs_theta"}};
    gStyle->SetPadRightMargin(0.15); 

    c1->DivideSquare(6);

    c1->cd(1);
    hKinMeasured->Draw("colz");

    c1->cd(2);
    hKinSampled->Draw("colz");

    c1->cd(3);
    auto eff {new TEfficiency(*hThetaCMRec, *hThetaCMSampled)};
    eff->Draw("AP");
    gPad->Update();

    c1->cd(4);
    hThetaCMSampled->Draw();


    c1->cd(5);
    hEx->Draw();

    c1->cd(6);
    hRange->Draw();


    c1->SaveAs(TString::Format("Graficas/Completo_Ex%.2f_incIdx%i.pdf", ExOg,incIdx));

    //eff->GetPaintedGraph()->GetXaxis()->SetTitle("#theta_{CM}");
    //eff->GetPaintedGraph()->GetYaxis()->SetTitle("Eficiencia");

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
    auto *c3{new TCanvas{"c3", "T3_vs_theta"}};
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
        c1->cd(4);
        hEx->Draw();

    c3->SaveAs(TString::Format("Graficas/Impacts/Impacts_Ex%.2f_incIdx%i.pdf", ExOg,incIdx));




}   
//####################################################### 