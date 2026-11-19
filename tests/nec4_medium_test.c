#include "nec4_medium.h"

#include <assert.h>
#include <complex.h>
#include <math.h>

static int close_complex(double complex a,double complex b,double tolerance)
{ return cabs(a-b)<=tolerance*fmax(1.,fmax(cabs(a),cabs(b))); }

int main(void)
{
  double x1[]={0.,0.},y1[]={0.,0.},z1[]={1.,0.};
  double x2[]={0.,0.},y2[]={0.,0.},z2[]={0.,-1.};
  nec4_medium_geometry geometry={2,x1,y1,z1,x2,y2,z2};
  nec4_medium_state state;
  double complex upper,lower,epsilon;
  nec4_medium_state_init(&state);
  assert(nec4_medium_classify(&geometry,4.,0.,100.,0,&state));
  assert(state.upper_segments==1&&state.lower_segments==1);
  assert(state.crossing_junctions==1&&state.interface_endpoints==2);
  assert(close_complex(state.refractive_index,2.,1e-14));
  assert(close_complex(state.intrinsic_impedance_ohm,376.730313668/2.,1e-14));
  assert(nec4_medium_endpoint_factor(&state,0,1,&upper));
  assert(nec4_medium_endpoint_factor(&state,1,0,&lower));
  assert(close_complex(upper/lower,1./4.,1e-14));
  assert(nec4_medium_update_frequency(&state,150.,9.-.3*I));
  assert(close_complex(state.epsilon_relative,9.-.3*I,1e-14));
  assert(nec4_medium_endpoint_factor(&state,0,1,&upper));
  assert(nec4_medium_endpoint_factor(&state,1,0,&lower));
  assert(close_complex(upper,1.,1e-14));
  assert(close_complex(lower,9.-.3*I,1e-14));
  assert(!nec4_medium_update_frequency(&state,150.,9.+.3*I));
  nec4_medium_state_free(&state);

  /* Endpoint classification is scale-aware and deterministic around z=0.
     Both interface endpoints below are inside the documented tolerance and
     are clustered into one crossing. */
  {
    double tx1[]={0.,0.},ty1[]={0.,0.},tz1[]={1.,0.};
    double tx2[]={0.,0.},ty2[]={0.,0.},tz2[]={1.e-14,-1.};
    nec4_medium_geometry tolerant={2,tx1,ty1,tz1,tx2,ty2,tz2};
    assert(nec4_medium_classify(&tolerant,4.,0.,100.,0,&state));
    assert(state.crossing_junctions==1&&state.interface_endpoints==2);
    assert(state.interface_tolerance_m>1.e-14);
    nec4_medium_state_free(&state);
  }

  z1[0]=1.;z2[0]=-1.;geometry.count=1;
  assert(!nec4_medium_classify(&geometry,4.,0.,100.,0,&state));
  assert(state.unsplit_crossings==1);
  nec4_medium_state_free(&state);

  z1[0]=z2[0]=0.;
  assert(!nec4_medium_classify(&geometry,4.,0.,100.,0,&state));
  assert(state.in_plane_segments==1);
  nec4_medium_state_free(&state);

  epsilon=nec4_medium_epsilon(13.,.02,10.,0);
  assert(cimag(epsilon)<0.&&close_complex(epsilon,
      13.-I*.02/(2.*acos(-1.)*10.e6*8.8541878128e-12),1e-14));
  assert(close_complex(nec4_medium_epsilon(4.,.25,10.,1),4.-.25*I,1e-14));
  nec4_medium_state_free(&state);
  return 0;
}
