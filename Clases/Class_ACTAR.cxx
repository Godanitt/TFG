#include "ActSRIM.h"
#include "ActKinematics.h"
#include "ActGeometry.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

#include "Rtypes.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
class ActarTPC
{
    double xActar{256}; 
    double yActar{256}; // Esta en mm
    double zActar{256};
    double zSil{1.5}; // Esta en mm
public:
    double get_xActar() { return xActar; };
    double get_yActar() { return yActar; }; 
    double get_zActar() { return zActar; };
    double get_zSil() { return zSil; };

};