/*
 * Optional M4 provenance harness. Build this against an unmodified DIRECTFN
 * checkout; it is not part of NEC2C and introduces no runtime dependency.
 */
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>

#include "directfn_algorithm_ea.h"
#include "directfn_algorithm_st.h"
#include "directfn_algorithm_va.h"
#include "directfn_contour.h"
#include "directfn_kernel_quad_scal.h"

using Directfn::Quadrilateral_EA;
using Directfn::Quadrilateral_ST;
using Directfn::Quadrilateral_VA;
using Directfn::QuadrilateralKernel_PlanarScalar;
using Directfn::SingularContour3xn;

template <typename Algorithm>
static void print_value(const char *name, SingularContour3xn &contour,
                        double k0, size_t order)
{
  std::unique_ptr<Algorithm> evaluator(new Algorithm());
  evaluator->set_wavenumber(k0);
  evaluator->set_Gaussian_orders_4(order, order, order, order);
  evaluator->set(contour);
  evaluator->calc_Iss();
  std::cout << name << ',' << std::setprecision(17)
            << evaluator->Iss()[0].real() << ',' << evaluator->Iss()[0].imag()
            << '\n';
}

int main()
{
  const double k0 = 2.1;
  const size_t order = 32;
  double r1[] = {0.0, 0.0, 0.0};
  double r2[] = {0.2, 0.0, 0.0};
  double r3[] = {0.2, 0.3, 0.0};
  double r4[] = {0.0, 0.3, 0.0};
  double r5_edge[] = {0.0, 0.6, 0.0};
  double r6_edge[] = {0.2, 0.6, 0.0};
  double r5_vertex[] = {0.4, 0.3, 0.0};
  double r6_vertex[] = {0.4, 0.6, 0.0};
  double r7_vertex[] = {0.2, 0.6, 0.0};
  SingularContour3xn self, edge, vertex;
  self.set_points(r1, r2, r3, r4);
  edge.set_points(r1, r2, r3, r4, r5_edge, r6_edge);
  vertex.set_points(r1, r2, r3, r4, r5_vertex, r6_vertex, r7_vertex);
  print_value<Quadrilateral_ST<QuadrilateralKernel_PlanarScalar>>(
      "same", self, k0, order);
  print_value<Quadrilateral_EA<QuadrilateralKernel_PlanarScalar>>(
      "common-edge", edge, k0, order);
  print_value<Quadrilateral_VA<QuadrilateralKernel_PlanarScalar>>(
      "common-vertex", vertex, k0, order);
  return 0;
}
