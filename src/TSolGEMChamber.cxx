#include <iostream>

#include "TSolGEMChamber.h"
#include "TSolGEMPlane.h"
#include "THaEvData.h"
#include "THaApparatus.h"
#include "TMath.h"
#include "ha_compiledata.h"

using namespace std;

TSolGEMChamber::TSolGEMChamber( const char *name, const char *desc )
  : THaDetector (name, desc)
{
  // For now at least we just hard wire two chambers
  fNPlanes = 2;
  fPlanes = new TSolGEMPlane*[fNPlanes];
  fChamberType = -1;
  return;
}

TSolGEMChamber::~TSolGEMChamber()
{
  for (UInt_t i = 0; i < fNPlanes; ++i)
    delete fPlanes[i];
  delete[] fPlanes;
  delete fShape;
}


const char* TSolGEMChamber::GetDBFileName() const {
    THaApparatus *app = GetApparatus();
    if( app )
      return Form ("%s.", app->GetName());
    else
      return fPrefix;
}

Int_t
TSolGEMChamber::ReadDatabase (const TDatime& date)
{
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  Int_t err = ReadGeometry (file, date, false);

  fclose(file);
  if (err)
    return err;

  err = InitPlane (0, TString (GetName()) + "x", TString (GetTitle()) +" x");
  if( err != kOK ) return err;
  err = InitPlane (1, TString (GetName()) + "y", TString (GetTitle()) +" y");
  if( err != kOK ) return err;

  return kOK;
}

Int_t
TSolGEMChamber::ReadGeometry (FILE* file, const TDatime& date,
			      Bool_t required)
{

  Int_t err = THaDetector::ReadGeometry (file, date, required);
  if (err)
    return err;

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

  fChamberType = type;
  
  if (fChamberType == 0){
    
    fShape = new TSolWedge;
    
    Double_t r0 = -999.0;
    Double_t r1 = -999.0;
    Double_t phi0 = -999.0;
    Double_t dphi = -999.0;
    Double_t z0 = -999.0;
    Double_t depth = -999.0;
    const DBRequest request[] =
      {
        {"r0",          &r0,           kDouble, 0, 1},
        {"r1",          &r1,           kDouble, 0, 1},
        {"phi0",        &phi0,         kDouble, 0, 1},
        {"dphi",        &dphi,         kDouble, 0, 1},
        {"z0",          &z0,           kDouble, 0, 1},
        {"depth",       &depth,        kDouble, 0, 1},
        {0}
      };
    err = LoadDB (file, date, request, fPrefix);

    if (err)
      return err;
  
    // Database specifies angles in degrees, convert to radians
    Double_t torad = atan(1) / 45.0;
    phi0 *= torad;
    dphi *= torad;

    fShape->SetGeometry (r0, r1, phi0, dphi);
    
    fOrigin[0] = (fShape->GetOrigin())[0];
    fOrigin[1] = (fShape->GetOrigin())[1];
    fOrigin[2] = z0;
    fSize[0] = (fShape->GetSize())[0];
    fSize[1] = (fShape->GetSize())[1];
    fSize[2] = depth;

    return kOK;
  }else if (fChamberType == 1) {
    fShape = new TSolRect;
    
    Double_t x0 = -999.0;
    Double_t x1 = -999.0;
    Double_t y0 = -999.0;
    Double_t y1 = -999.0;
    Double_t z0 = -999.0;
    Double_t depth = -999.0;
    
    const DBRequest request[] =
      {
        {"x0",          &x0,           kDouble, 0, 1},
        {"x1",          &x1,           kDouble, 0, 1},
        {"y0",          &y0,           kDouble, 0, 1},
        {"y1",          &y1,           kDouble, 0, 1},
        {"z0",          &z0,           kDouble, 0, 1},
        {"depth",       &depth,        kDouble, 0, 1},
        {0}
      };
    err = LoadDB (file, date, request, fPrefix);

    if (err)
      return err;
      
    fShape->SetGeometry(x0, x1, y0, y1);
    
    fOrigin[0] = (fShape->GetOrigin())[0];
    fOrigin[1] = (fShape->GetOrigin())[1];
    fOrigin[2] = z0;
    fSize[0] = (fShape->GetSize())[0];
    fSize[1] = (fShape->GetSize())[1];
    fSize[2] = depth;
    
  }else{
    return -1;
  }
  
  return kOK;
}


Int_t
TSolGEMChamber::Decode (const THaEvData& ed )
{
  for (UInt_t i = 0; i < GetNPlanes(); ++i)
    {
      GetPlane (i).Decode (ed);
    }
  return 0;
}

Int_t
TSolGEMChamber::InitPlane (const UInt_t i, const char* name, const char* desc)
{
  
  fPlanes[i] = new TSolGEMPlane (name, desc, this);
  fPlanes[i]->SetName (name);
  return fPlanes[i]->Init();
}

void
TSolGEMChamber::Print (const Bool_t printplanes)
{
  cout << "I'm a GEM chamber named " << GetName() << endl;
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

  fShape->Print();
  /*if (fChamberType == 0){
  cout << "  Wedge geometry: r0: " << fShape->GetR0()
       << " r1: " << fShape->GetR1()
       << " phi0: " << fShape->GetPhi0()*TMath::RadToDeg()
       << " dphi: " << fShape->GetDPhi()*TMath::RadToDeg()
       << endl;
  }else if (fChamberType == 1){
    cout << "  Rectangle geometry: x0: " << fShape->GetR0()
       << " x1: " << fShape->GetR1()
       << " y0: " << fShape->GetPhi0()
       << " y1: " << fShape->GetDPhi()
       << endl;
  }*/
  if (printplanes)
    for (UInt_t i = 0; i < GetNPlanes(); ++i)
      {
	fPlanes[i]->Print();
      }
}
