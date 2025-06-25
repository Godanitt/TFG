
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
void Sampleo() {
    
gStyle->SetLabelSize(0.05,"XYZ");  // Tamaño de los números en el eje X
gStyle->SetLabelFont(62,"XYZ");
gStyle->SetTitleSize(0.05,"XYZ");  // Tamaño de los números en el eje X
gStyle->SetTitleFont(62);      // Sin segundo argumento → título general
gStyle->SetTitleFont(62,"XYZ");
gStyle->SetTitleSize(0.05); 
// Margenes: 
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadRightMargin(0.02);
gStyle->SetPadTopMargin(0.05);
gStyle->SetPadBottomMargin(0.14);

gROOT->ForceStyle();


auto* kinEx0 {new ActPhysics::Kinematics("11Li", "d", "t",11 * 7.5, 0.0)}; // cinemática
auto* kinEx2 {new ActPhysics::Kinematics("11Li", "d", "t",11 * 7.5, 0.2)}; // cinemática


auto* g3Li10Ex0 {kinEx0->GetKinematicLine3()};
auto* g4Li10Ex0 {kinEx0->GetKinematicLine4()};

auto* g3Li10Ex2 {kinEx2->GetKinematicLine3()};
auto* g4Li10Ex2 {kinEx2->GetKinematicLine4()};



auto* canvas10 {new TCanvas {"canvas10", "Li10",800,400}};


// Grosor de los ejes
gStyle->SetLineWidth(2);      // Más grueso el eje X
gStyle->SetLineWidth(2);      // Más grueso el eje Y

canvas10->DivideSquare(2);
canvas10->cd(1);

g3Li10Ex0->SetLineColor(kBlue);
g3Li10Ex0->SetLineWidth(1.5);
g3Li10Ex2->SetLineColor(kRed);
g3Li10Ex2->SetLineWidth(1.5);

g3Li10Ex0->Draw("al");
g3Li10Ex2->Draw("l same");

auto legend = new TLegend(0.65,0.8,0.98,0.95);
legend->SetTextSize(0.05);  // Tamaño del texto
legend->SetTextFont(62);     // Fuente en negrita



legend->AddEntry(g3Li10Ex0,"E_{x}=0.0 MeV","l");
legend->AddEntry(g3Li10Ex2,"E_{x}=0.2 MeV","l");
legend->Draw();
g3Li10Ex0->SetTitle(";#theta_{3Lab} [#circ];E_{3Lab} [MeV]");

canvas10->RedrawAxis();         // Redibuja los ejes, opcional
canvas10->cd(2);


canvas10->Update();
canvas10->SaveAs(TString::Format("Graficas/Cinematica.pdf"));

canvas10->cd(2);

g4Li10Ex0->SetLineColor(kBlue);
g4Li10Ex0->SetLineWidth(2);
g4Li10Ex2->SetLineColor(kRed);
g4Li10Ex2->SetLineWidth(2);

g4Li10Ex0->Draw("al");
g4Li10Ex2->Draw("l same");

auto legend2= new TLegend(0.65,0.8,0.98,0.95);
legend2->SetTextSize(0.05);  // Tamaño del texto
legend2->SetTextFont(62);     // Fuente en negrita

legend2->AddEntry(g4Li10Ex0,"E_{x}=0.0 MeV","l");
legend2->AddEntry(g4Li10Ex2,"E_{x}=0.2 Mev","l");
legend2->Draw();
g4Li10Ex0->SetTitle(TString::Format(";#theta_{4Lab} [#circ];E_{4Lab} [MeV]"));

canvas10->RedrawAxis();         // Redibuja los ejes, opcional

canvas10->RedrawAxis();
canvas10->cd(2);


canvas10->Update();
canvas10->SaveAs(TString::Format("/home/daniel/GitHub/TFG/Memoria/Imagenes/Cinematica.pdf"));

}