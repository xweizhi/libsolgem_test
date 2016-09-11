TSolSpec* dds;
TSolGEMChamber* ddy;
TSolSimGEMDigitization* ddd;
TSolGEMData* ddgd;

void DigDemo()
{
  gSystem->Load("../libsolgem.so");

  dds = new TSolSpec ("spectrometer", "SOLiD spectrometer");
  dds->Init();

  ddy = new TSolGEMChamber ("testchamber","Test chamber");
  ddy->SetName("testchamber");
  ddy->Init();
  ddy->Print();

  dds->AddGEM (ddy);

  ddd = new TSolSimGEMDigitization (*dds);
  
  ddgd = new TSolGEMData (1); // 1 hit wonder
  ddgd->SetRun (1000);
  ddgd->SetEvent (0);
  ddgd->SetNHit(1);
  ddgd->SetMomentum (0, TVector3 (20, 20, 3000.00));
  ddgd->SetHitEntrance (0, TVector3 (0.27, 0.00, 1.55) * 1000.0); // mm
  ddgd->SetHitExit (0, TVector3 (0.2809, 0.000, 1.58) * 1000.0);
  ddgd->SetHitReadout (0, TVector3 (0.2846, 0.000, 1.59) * 1000.0);
  ddgd->SetHitTime(0, 0.0); //ns
  ddgd->SetHitEnergy (0, 100); // eV
  ddgd->SetHitChamber (0, 0);
  ddgd->SetParticleID (0, 1);
  ddgd->SetParticleType (0, 0);
  //ddgd->SetNHit(1);
  ddgd->Print();
  ddgd->PrintHit (0);

  ddd->InitTree (*dds, "digdemo.root");
  ddd->Digitize (*ddgd, *dds);
  ddd->PrintCharges();
  ddd->WriteTree();
  ddd->CloseTree();
}
