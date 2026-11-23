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
  assert(state.one_sided_interface_ends==0);
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

  /* The published radial-screen topology has one upper monopole ray and
     several lower-medium radials meeting at the interface. */
  {
    double rx1[]={0.,0.,0.,0.,0.,0.};
    double ry1[]={0.,0.,0.,0.,0.,0.};
    double rz1[]={1.,0.,0.,0.,0.,0.};
    double rx2[]={0.,1.,.309016994,-.809016994,-.809016994,.309016994};
    double ry2[]={0.,0.,.951056516,.587785252,-.587785252,-.951056516};
    double rz2[]={0.,-1.e-5,-1.e-5,-1.e-5,-1.e-5,-1.e-5};
    nec4_medium_geometry radial={6,rx1,ry1,rz1,rx2,ry2,rz2};
    assert(nec4_medium_classify(&radial,10.,1.,100.,1,&state));
    assert(state.crossing_junctions==1&&state.radial_screen_junctions==1&&
      state.upper_segments==1&&state.lower_segments==5&&
      state.unsupported_interface_junctions==0);
    assert(nec4_medium_endpoint_factor(&state,0,1,&upper));
    for(int radial_index=1;radial_index<6;++radial_index)
    {
      assert(nec4_medium_endpoint_factor(&state,(size_t)radial_index,0,&lower));
      assert(close_complex(upper/lower,1./(10.-I),1e-14));
    }
    nec4_medium_state_free(&state);
  }

  /* A generic branched interface is not silently treated as a radial
     screen: the lower rays must form a balanced, nearly horizontal fan. */
  {
    double bx1[]={0.,0.,0.,0.},by1[]={0.,0.,0.,0.};
    double bz1[]={1.,0.,0.,0.};
    double bx2[]={0.,1.,0.,-0.5},by2[]={0.,0.,1.,0.};
    double bz2[]={0.,-1.e-5,-1.e-5,-1.e-5};
    nec4_medium_geometry branch={4,bx1,by1,bz1,bx2,by2,bz2};
    assert(!nec4_medium_classify(&branch,10.,1.,100.,1,&state));
    assert(state.radial_screen_junctions==0&&
      state.unsupported_interface_junctions==1);
    nec4_medium_state_free(&state);
  }


  /* A single free end on the interface is a supported non-contacting
     topology, distinct from a crossing junction or an in-plane segment. */
  {
    double fx1[]={0.},fy1[]={0.},fz1[]={0.};
    double fx2[]={0.},fy2[]={0.},fz2[]={1.};
    nec4_medium_geometry free_end={1,fx1,fy1,fz1,fx2,fy2,fz2};
    assert(nec4_medium_classify(&free_end,4.,0.,100.,0,&state));
    assert(state.upper_segments==1&&state.interface_endpoints==1&&
      state.one_sided_interface_ends==1&&state.crossing_junctions==0&&
      state.unsupported_interface_junctions==0);
    nec4_medium_state_free(&state);
  }

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
