#include "TSolGEMPlane.h"

#include "TClonesArray.h"

#include "TSolGEMChamber.h"
#include "TSolGEMCluster.h"
#include "TSolWedge.h"
#include "THaEvData.h"
#include "TMath.h"
#include "ha_compiledata.h"

#include <iostream>
#include <cassert>

using namespace std;

TSolGEMPlane::TSolGEMPlane()
  : THaSubDetector(), fShape(0)
{
  //  fClusters = new TClonesArray("TSolGEMCluster", 100);  
  return;
}

TSolGEMPlane::TSolGEMPlane( const char *name, const char *desc,
			    THaDetectorBase* parent )
  : THaSubDetector (name, desc, parent)
{
  //  fClusters = new TClonesArray("TSolGEMCluster", 100);  
  return;
}

TSolGEMPlane::~TSolGEMPlane()
{
  //  delete fClusters;
  //delete fShape;
}

Int_t 
TSolGEMPlane::ReadDatabase (const TDatime& date)
{
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  Int_t err = ReadGeometry (file, date, false);

  fclose(file);
  if (err)
    return err;

  return kOK;
}

Int_t 
TSolGEMPlane::ReadGeometry (FILE* file, const TDatime& date,
			    Bool_t required)
{
  // Get x/y position, size, and angles from database if and only
  // if parent is null otherwise copy from parent

  // Note that origin is in lab frame, size is in wedge frame.

  Int_t err;
  Double_t torad = atan(1) / 45.0;

  TSolGEMChamber* parent = (TSolGEMChamber*) GetParent();
  Double_t z0;
  Double_t depth;
  
  if (parent != NULL)
    {
      fOrigin = parent->GetOrigin();
      fSize[0] = (parent->GetSize())[0];
      fSize[1] = (parent->GetSize())[1];
      fSize[2] = (parent->GetSize())[2];
      
      if (parent->GetChamberType() == 0) {
        fShape = new TSolWedge;
      }else if (parent->GetChamberType() == 1){
        fShape = new TSolRect;
      }
         
     fShape->SetGeometry (parent->GetShape());
     
      z0 = fOrigin[2];
      depth = fSize[2];
    }
  else
    {
      
      
      Int_t type = -1;
  const DBRequest typeRequest[] = 
  {
    {"type",         &type,         kInt,    0,  1},
    {0}
  };
  

  err = LoadDB (file, date, typeRequest, fPrefix);
  if (err) return err;
  
  if (!(type == 0 || type == 1)) {
    cout<<"unknown GEM geometric type "<<type<<endl;
    return -1;
  }
  
  
     
  if (type == 0){
    fShape = new TSolWedge;
  
      Double_t r0 = -999.0;
      Double_t r1 = -999.0;
      Double_t phi0 = -999.0;
      Double_t dphi = -999.0;
  
      const DBRequest request[] = 
	{
	  {"r0",          &r0,           kDouble, 0, 1},
	  {"r1",          &r1,           kDouble, 0, 1},
	  {"phi0",        &phi0,         kDouble, 0, 1},
	  {"dphi",        &dphi,         kDouble, 0, 1},
	  {0}
	};
      err = LoadDB( file, date, request, fPrefix );
      
      if (err)
	return err;
      // Database specifies angles in degrees, convert to radians
      phi0 *= torad;
      dphi *= torad;
      
      fShape->SetGeometry (r0, r1, phi0, dphi);

      fOrigin[0] = (fShape->GetOrigin())[0];
      fOrigin[1] = (fShape->GetOrigin())[1];
      fSize[0] = (fShape->GetSize())[0];
      fSize[1] = (fShape->GetSize())[1];
      z0 = -999.0;
      depth = -999.0;
    }
  else if (type == 1) {
    fShape = new TSolRect;
    
    Double_t x0 = -999.0;
    Double_t x1 = -999.0;
    Double_t y0 = -999.0;
    Double_t y1 = -999.0;
    
    
    const DBRequest request[] =
      {
        {"x0",          &x0,           kDouble, 0, 1},
        {"x1",          &x1,           kDouble, 0, 1},
        {"y0",          &y0,           kDouble, 0, 1},
        {"y1",          &y1,           kDouble, 0, 1},
        {0}
      };
    err = LoadDB (file, date, request, fPrefix);

    if (err)
      return err;
      
    fShape->SetGeometry(x0, x1, y0, y1);
    
    fOrigin[0] = (fShape->GetOrigin())[0];
    fOrigin[1] = (fShape->GetOrigin())[1];
    
    fSize[0] = (fShape->GetSize())[0];
    fSize[1] = (fShape->GetSize())[1];
    
  }else{
    return -1;
  }
  }
  
  
  
  const DBRequest request[] = 
    {
      {"stripangle",  &fSAngle,      kDouble, 0, 1},
      {"pitch",       &fSPitch,      kDouble, 0, 1},
      {"z0",          &z0,           kDouble, 0, 1},
      {"depth",       &depth,        kDouble, 0, 1},
      {0}
    };
  err = LoadDB( file, date, request, fPrefix );
  
  if (err)
    return err;

  fSAngle *= torad;

  SetRotations();
  fOrigin[2] = z0;
  fSize[2] = depth;
  
  // Get numbers of strips

  Double_t xs0 = (GetSize())[0];
  Double_t ys0 = (GetSize())[1];
  Double_t xs[4] = {xs0, xs0, -xs0, -xs0};
  Double_t ys[4] = {ys0, -ys0, ys0, -ys0};
  Double_t smin = 1e9;
  Double_t smax = 1e-9;
  for (UInt_t i = 0; i < 4; ++i)
    {
      PlaneToStrip (xs[i], ys[i]);
      Double_t s = (xs[i] / GetSPitch());
      if (s < smin) smin = s;
      if (s > smax) smax = s;
    }
  fNStrips = (Int_t)(smax - smin);
  fSBeg = -fNStrips * fSPitch * 0.5;
  cout<<fNStrips<<" "<<fSBeg<<endl;
  return kOK;
}

Int_t TSolGEMPlane::Decode( const THaEvData &d ){
    // Clusters get made as so

  //    int i = 0;

    //    new ((*fClusters)[i]) TSolGEMCluster();

    return 0;
}

Double_t 
TSolGEMPlane::GetSAngle()   const
{
  return fSAngle;
}

void
TSolGEMPlane::LabToStrip (Double_t& x, Double_t& y) const
{
  x -= (GetOrigin())[0];
  y -= (GetOrigin())[1];
  Double_t temp = x;
  x = fCLS * x - fSLS * y;
  y = fSLS * temp + fCLS * y;
  return;
}

void
TSolGEMPlane::StripToPlane (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCWS * x + fSWS * y;
  y = -fSWS * temp + fCWS * y;
  return;
}

void
TSolGEMPlane::StripToLab (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCLS * x + fSLS * y;
  y = -fSLS * temp + fCLS * y;
  x += (GetOrigin())[0];
  y += (GetOrigin())[1];
  return;
}

Double_t TSolGEMPlane::StripNumtoStrip( Int_t strip )
{
    // Gives x coordinate in strip frame of a wire
    return (strip - GetStrip(0.,0.))*GetSPitch();
}


Double_t TSolGEMPlane::StriptoProj( Double_t s )
{
    //this is probably wrong? (r1+r2)/2 is not the same as the center of the bounding box
    
    // Gives coordinate in projection frame from strip frame x
    Double_t r = sqrt(pow((GetOrigin())[0],2) + pow((GetOrigin())[1],2));
    return s + r*fCWS;
}


Double_t TSolGEMPlane::StripNumtoProj( Int_t s ){
    // Gives coordinate in projection frame from strip number
    return StriptoProj( StripNumtoStrip(s) );
}

Double_t 
TSolGEMPlane::GetStripLowerEdge (UInt_t is) const {return (fSBeg + is * GetSPitch());}

Double_t 
TSolGEMPlane::GetStripUpperEdge (UInt_t is) const {return GetStripLowerEdge (is) + GetSPitch();}


Int_t
TSolGEMPlane::GetStripUnchecked( Double_t x ) const
{
  // Get strip number for given x-coordinate in strip frame,
  // no questions asked. Caller must check return value

  return (Int_t) ((x - fSBeg) / GetSPitch());
}

Int_t
TSolGEMPlane::GetStripInRange( Double_t x ) const
{
  // Get strip number for given x-coordinate in strip frame
  // and, if out of range, limit it to allowable values.

  Int_t s = GetStripUnchecked(x);
  if( s < 0 )              s = 0;
  if( s >= GetNStrips() )  s = GetNStrips()-1;
  return s;
}
    
Int_t
TSolGEMPlane::GetStrip (Double_t x, Double_t yc) const
{
  // Strip number corresponding to coordinates x, y in 
  // strip frame, or -1 if outside (2-d) bounds

  Double_t xc = x;
  StripToLab (xc, yc);

  if (!fShape->Contains (xc, yc))
    return -1;

  Int_t s = GetStripInRange(x);
  assert( s >= 0 && s < GetNStrips() ); // by construction in ReadGeometry()
  return s;
}

void 
TSolGEMPlane::Print() const
{
  cout << "I'm a GEM plane named " << GetName() << endl;

  TVector3 o (GetOrigin());
  cout << "  Origin: " << o(0) << " " << o(1) << " " << o(2)
       << " (rho,theta,phi)=(" << o.Mag() << "," << o.Theta()*TMath::RadToDeg()
       << "," << o.Phi()*TMath::RadToDeg() << ")"
       << endl;

#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
  const Double_t* s = GetSize();
#else
  const Float_t* s = GetSize();
#endif
  cout << "  Size:   " << s[0] << " " << s[1] << " " << s[2] << endl;
   
  dynamic_cast<TSolGEMChamber*>(GetParent())->GetShape()->Print();   
 
  cout << "  " << GetNStrips() << " strips"
       << ", angle " << GetSAngle()*TMath::RadToDeg()
       << ", start " << fSBeg << " " << 0.5*(GetStripLowerEdge(0)+GetStripUpperEdge(0))
       << ", pitch " << GetSPitch()
       << endl;
}

void
TSolGEMPlane::SetRotations()
{
  // Set rotation angle trig functions
  //cout<<"strip angles: "<<GetSAngle()<<" "<<GetAngle()<<endl;
  fCWS = cos (-GetSAngle());
  fCLS = cos (-GetAngle()-GetSAngle());
  fSWS = sin (-GetSAngle());
  fSLS = sin (-GetAngle()-GetSAngle());
}
