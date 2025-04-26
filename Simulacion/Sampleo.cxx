
# include "../Clases/Class_Colision.cxx"
# include "../Clases/Estilo_Publicacion.cxx"
# include "../Clases/Class_ACTAR.cxx"
/*

En esta vamos a crear una función que nos va a dibuar las gráficas de T3,T4 vs theta3, theta4.
Entrada: 
- Particula p1: la partícula incidente
- Particula p2: la particula blanco
- Particula p3: la partícula producto ligera
- Particula p4: la partícula producto pesada
- ActPhysics::Kinematics Li10: el objeto de cinemática que nos da la información de la colisión
- double tbeam: energía del haz (en MeV) (por nucleon)

Posible implementaciones para el futuro:
- Nombre pdf
*/   
void SampleoCinematica(Particula p1, Particula p2, Particula  p3, Particula p4, ActPhysics::Kinematics* Li10, double tbeam, double Ex) {

auto* gother10 {new TGraph};
auto* gother3 {new TGraph};
auto* g3Li10 {Li10->GetKinematicLine3()};
auto* g4Li10 {Li10->GetKinematicLine4()};



// Queremos representar todas las posibilidades -> Recorremos todos los ángulos
for(double thetaCM = 0; thetaCM <= 180; thetaCM += 1)
{       
    auto [T3,Theta3] {valores3(p1, p2, p3, p4, p1.get_A()*tbeam,thetaCM * TMath::DegToRad(),Ex)};
    auto [T4,Theta4] {valores4(p1, p2, p3, p4, p1.get_A()*tbeam,thetaCM * TMath::DegToRad(),Ex)};
    // std::cout << "Theta4: " << Theta4 << " T4 : " << T4 << '\n';
    if(std::isfinite(Theta4) && std::isfinite(T4))
        gother10->SetPoint(gother10->GetN(), Theta4 * TMath::RadToDeg(), T4);
    if(std::isfinite(Theta3) && std::isfinite(T3))
        gother3->SetPoint(gother3->GetN(), Theta3 * TMath::RadToDeg(), T3);
}

auto* canvas10 {new TCanvas {"canvas10", "Li10"}};

SetEstiloPublicacion();
canvas10->DivideSquare(2);
canvas10->cd(1);
SetEstiloPublicacion();
g3Li10->SetLineColor(kBlue);
g3Li10->SetLineWidth(4);
gother3->SetLineColor(kGray+1);
gother3->SetLineWidth(1);
gother3->SetMarkerStyle(45);
gother3->SetMarkerColor(kRed);
gother3->SetMarkerSize(0.5);  // Un pelín más grandes
g3Li10->Draw("al");
gother3->Draw("l same"); 

auto legend = new TLegend(0.7,0.7,0.95,0.95);
legend->AddEntry(gother3,"Cinematica","l");
legend->AddEntry(g3Li10,"ActPhysics","l");
legend->Draw();
g3Li10->SetTitle("Cinematica 11Li(d,t);#theta_{3Lab} [#circ];E_{3Lab}");

canvas10->RedrawAxis();         // Redibuja los ejes, opcional

canvas10->RedrawAxis();
canvas10->cd(2);

SetEstiloPublicacion();
canvas10->Update();
canvas10->SaveAs(TString::Format("Graficas/Graficas T vs theta (sampleada)/Cinematica.pdf"));

canvas10->cd(2);
SetEstiloPublicacion();
g4Li10->SetLineColor(kBlue);
g4Li10->SetLineWidth(4);
gother10->SetLineColor(kGray+1);
gother10->SetLineWidth(1);
gother10->SetMarkerStyle(45);
gother10->SetMarkerColor(kRed);
gother10->SetMarkerSize(0.5);  // Un pelín más grandes
g4Li10->Draw("al");
gother10->Draw("l same"); 

auto legend2 = new TLegend(0.7,0.7,0.95,0.95);
legend2->AddEntry(gother10,"Cinematica","l");
legend2->AddEntry(g4Li10,"ActPhysics","l");
legend2->Draw();
g4Li10->SetTitle(TString::Format("Cinematica Sampleada Ex=%.2f;#theta_{3Lab} [#circ];E_{3Lab}",Ex));

canvas10->RedrawAxis();         // Redibuja los ejes, opcional

canvas10->RedrawAxis();
canvas10->cd(2);

SetEstiloPublicacion();
canvas10->Update();
canvas10->SaveAs(TString::Format("Graficas/Cinematica Sampleada/Cinematica_%.2fEx.pdf",Ex));

}