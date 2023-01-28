/*
The conversion from unitness to mV is 120 * Y0 - 81.0
*/
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{

 realtype Y[NEQ];
 realtype dY[NEQ];

 realtype IV;
 UserData data;
 data = (UserData) user_data;

 int i;

 IV = data->IV;

  for(i=0;i<NEQ;i++)
    Y[i] = Ith(y,i+1);

// FK basal parameters.
/* data->p[0] */ realtype tau_d      = data->p[0];  //   0.25;
/* data->p[1] */ realtype tau_r      = data->p[1];  //   33  ;
/* data->p[2] */ realtype tau_si     = data->p[2];  //   30  ;
/* data->p[3] */ realtype tau_0      = data->p[3];  //   12.5;
/* data->p[4] */ realtype tau_vp     = data->p[4];  //   3.33;
/* data->p[5] */ realtype tau_vm1    = data->p[5];  //   1250;
/* data->p[6] */ realtype tau_vm2    = data->p[6];  //   19.6;
/* data->p[7] */ realtype tau_wp     = data->p[7];  //   870 ;
/* data->p[8] */ realtype tau_wm     = data->p[8];  //   41  ;
/* data->p[9] */ realtype u_c        = data->p[9];  //   0.13;
/* data->p[10] */ realtype u_v       = data->p[10]; //   0.04;
/* data->p[11] */ realtype usi_c     = data->p[11]; //   0.85;
/* data->p[12] */ realtype k         = data->p[12]; //   1.0 ;

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

  realtype Jfi, Jso, Jsi, tau_vm;

  if (Y[0]<u_c) {
    Jfi    =  0.0;
    Jso    =  2*Y[0]/(tau_0);
    tau_vm = (Y[0] > u_v) ? (tau_vm1) : (tau_vm2);
    dY[1]  = (1.0 - Y[1])/tau_vm;
    dY[2]  = (1.0 - Y[2])/tau_wm;
  } else {
    Jfi    = - Y[1]/(tau_d)*(1.-Y[0])*(Y[0]-u_c);
    Jso    =   2.0 /(tau_r);
    dY[1]  = - Y[1]/tau_vp;
    dY[2]  = - Y[2]/tau_wp;
  }
    Jsi    = - Y[2]/tau_si*(1.0 + tanh(k*(Y[0]-usi_c)));
  dY[0]    = - Jfi - Jso - Jsi + IV + data->delta_V; // the data->delta_V is the diffusion current comes from the other cell to this cell.

 for(i=0;i<NEQ;i++)
  Ith(ydot,i+1) = dY[i];

data->dvdt = dY[0];

 return(0);
}

