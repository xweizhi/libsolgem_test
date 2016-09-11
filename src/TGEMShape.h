#ifndef TGEMSHAPE_H
#define TGEMSHAPE_H
#include <iostream>
#include <Rtypes.h>
#include <TMath.h>
class TGEMShape
{
  public:
  virtual ~TGEMShape();
  virtual void SetGeometry (TGEMShape* theShape) = 0;
  virtual void SetGeometry (Double_t a, Double_t b, Double_t c, Double_t d) = 0;
  virtual Bool_t Contains (Double_t x, Double_t y) const = 0;
  virtual void LabToWedge (Double_t& x, Double_t& y) const = 0;
  virtual void WedgeToLab (Double_t& x, Double_t& y) const = 0;
  virtual std::vector < Double_t > GetOrigin() const = 0;
  virtual std::vector < Double_t > GetSize() const = 0;
  virtual Double_t GetAngle() const = 0;
  virtual void Print () = 0; 
  protected:
  std::vector < Double_t > fOrigin;   // x, y
  std::vector < Double_t > fSize;     // x, y
  
};



#endif
