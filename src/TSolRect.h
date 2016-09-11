
#ifndef __TSOLRECT_H
#define __TSOLRECT_H

#include <cmath>
#include <vector>
#include "TGEMShape.h"
#include <Rtypes.h>

class TDatime;


class TSolRect : public TGEMShape
{
 public:
  TSolRect (Double_t xlow = -0.2752, Double_t xup = 0.2752, Double_t ylow = -0.6144, Double_t yup = 0.6144);
  virtual ~TSolRect() {};

  Double_t GetX1() const {return fX1;};
  Double_t GetX2() const {return fX2;};
  Double_t GetY1() const {return fY1;};
  Double_t GetY2() const {return fY2;};
  std::vector < Double_t > GetOrigin() const {return fOrigin;};
  std::vector < Double_t > GetSize() const {return fSize;};

  std::vector < Double_t > Bounds() const;
  Bool_t Contains (Double_t x, Double_t y) const;

  void SetGeometry (TGEMShape* theShape);
  void SetGeometry (const Double_t x1,
                    const Double_t x2,
                    const Double_t y1,
                    const Double_t y2);
    
  Double_t GetAngle() const {return fPhi; }; // rotation angle between lab and wedge frame

  // Frame conversions
  void LabToWedge (Double_t& x, Double_t& y) const;  // input and output in meters
  void WedgeToLab (Double_t& x, Double_t& y) const;  // input and output in meters
  void Print();
 private:
  void SetRotations();

  Double_t fX1;
  Double_t fX2;
  Double_t fY1;
  Double_t fY2;
  Double_t fPhi;

  //std::vector < Double_t > fOrigin;   // x, y
  //std::vector < Double_t > fSize;     // x, y

  // Trig functions for rotations
  Double_t fCLW; // cos (lab to wedge angle)
  Double_t fSLW; // sin...
};

#endif
