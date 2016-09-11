#include "TSolRect.h"

#include <iostream>
using namespace std;

#include <TDatime.h>

TSolRect::TSolRect (Double_t x1, Double_t x2, Double_t y1, Double_t y2)
  : fX1(-0.2752), fX2(0.2752), fY1(-0.6144), fY2(0.6144), fPhi(0.)
{
  SetGeometry (x1, x2, y1, y2);
}

std::vector < Double_t > 
TSolRect::Bounds() const
{
  // Return bounds: xmin, ymin, xmax, ymax
  // Each is either +-r1 or a coordinate of a corner.

  std::vector < Double_t > v(4); 

  v[0] = fX1;
  v[1] = fY1;
  v[2] = fX2;
  v[3] = fY2;

  return v;
}

Bool_t
TSolRect::Contains (Double_t x, Double_t y) const
{
  if (x <= fX2 && x >= fX1 && y <= fY2 && y >= fY1) return true;
  else return false;

}
//__________________________________________________________________
void TSolRect::SetGeometry(TGEMShape* theShape)
{
  Double_t x1 = dynamic_cast<TSolRect*>(theShape)->GetX1();
  Double_t x2 = dynamic_cast<TSolRect*>(theShape)->GetX2();
  Double_t y1 = dynamic_cast<TSolRect*>(theShape)->GetY1();
  Double_t y2 = dynamic_cast<TSolRect*>(theShape)->GetY2();
  SetGeometry(x1, x2, y1, y2);
}
//__________________________________________________________________
void
TSolRect::SetGeometry (const Double_t x1,
                        const Double_t x2,
                        const Double_t y1,
                        const Double_t y2)
{
  // Note that origin is in lab frame, size is in wedge frame.
  fX1 = x1;
  fX2 = x2;

  fY1 = y1;
  fY2 = y2;
  vector <Double_t> bounds = Bounds();

  SetRotations();

  // Origin is center of bounding box, converted to lab frame
  fOrigin.resize (2);
  
  Double_t xc = (bounds[2]+bounds[0])/2;
  Double_t yc = (bounds[3]+bounds[1])/2;
  
  fOrigin[0] = xc;
  fOrigin[1] = yc;
  
  // Size is half size of bounding box, in wedge frame
  fSize.resize (2);
  fSize[0] = (bounds[2] - bounds[0]) / 2.;
  fSize[1] = (bounds[3] - bounds[1]) / 2.;
}

void
TSolRect::LabToWedge (Double_t& x, Double_t& y) const
{
  x -= (GetOrigin())[0];
  y -= (GetOrigin())[1];
  Double_t temp = x;
  x = fCLW * x - fSLW * y;
  y = fSLW * temp + fCLW * y;
  return;
}

void
TSolRect::WedgeToLab (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCLW * x + fSLW * y;
  y = -fSLW * temp + fCLW * y;
  x += (GetOrigin())[0];
  y += (GetOrigin())[1];
  return;
}

void
TSolRect::SetRotations()
{
  // Set rotation angle trig functions

  fCLW = cos (-GetAngle());
  fSLW = sin (-GetAngle());
}

//______________________________________________________
void TSolRect::Print()
{
   cout << "  Rectangle geometry: x0: " << GetX1()
       << " x1: " << GetX2()
       << " y0: " << GetY1()
       << " y1: " << GetY2()
       << endl;

}

