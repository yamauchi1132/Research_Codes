#pragma once

class Profile
{
public:
  clock_t start, end, log1, log2;
  double Init, Init_dens, Init_hydro, InitialKick, FullDrift, AdjustPosition, Predict, Decompose, Exchange, DensCalc, SetPres, HydroForceCalc, FinalKick, Check, All, WriteFile; 
 
  double time_calc(clock_t start, clock_t end)
  {
    double buf = (double)(end - start)/CLOCKS_PER_SEC;
    return buf;
  }

  void disp()
  {
    fprintf(stdout, "\n######### Time Profile #########\n");
    fprintf(stdout, "All time = %lf second\n", All);
    fprintf(stdout, "\n");
    fprintf(stdout, "Init_dens time = %lf second\n", Init_dens);
    fprintf(stdout, "Init_hydro time = %lf second\n", Init_hydro);
    fprintf(stdout, "Init All time = %lf second\n", Init);
    fprintf(stdout, "\n");
    fprintf(stdout, "InitialKick time = %lf second\n", InitialKick);
    fprintf(stdout, "FullDrift time = %lf second\n", FullDrift);
    fprintf(stdout, "AdjustPosition time = %lf second\n", AdjustPosition);
    fprintf(stdout, "Predict time = %lf second\n", Predict);
    fprintf(stdout, "Decompose time = %lf second\n", Decompose);
    fprintf(stdout, "Exchange time = %lf second\n", Exchange);
    fprintf(stdout, "DensCalc time = %lf second\n", DensCalc);
    fprintf(stdout, "SetPres time = %lf second\n", SetPres);
    fprintf(stdout, "HydroForceCalc time = %lf second\n", HydroForceCalc);
    fprintf(stdout, "FinalKick time = %lf second\n", FinalKick);
    fprintf(stdout, "WriteFile time = %lf second\n", WriteFile);
    fprintf(stdout, "Check time = %lf second\n", Check);
    fprintf(stdout, "#################################\n\n");
  }

  Profile(){
    InitialKick = 0;
    FullDrift = 0;
    AdjustPosition = 0;
    Predict =0;
    Decompose = 0;
    Exchange = 0;
    DensCalc = 0;
    SetPres = 0;
    HydroForceCalc = 0;
    FinalKick = 0;
    Check = 0;
    WriteFile = 0;
  }
};
