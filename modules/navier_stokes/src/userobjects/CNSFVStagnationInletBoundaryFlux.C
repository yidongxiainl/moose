/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CNSFVStagnationInletBoundaryFlux.h"

template <>
InputParameters
validParams<CNSFVStagnationInletBoundaryFlux>()
{
  InputParameters params = validParams<BoundaryFluxBase>();

  params.addClassDescription(
    "A user objec that computes the stagnation inlet boundary flux.");

  params.addRequiredParam<UserObjectName>(
    "bc_uo",
    "Name for boundary condition user object");

  params.addRequiredParam<UserObjectName>(
    "fluid_properties",
    "Name for fluid properties user object");

  params.addRequiredParam<Real>(
    "stagnation_inlet_pressure",
    "Stagnation inlet pressure");

  params.addRequiredParam<Real>(
    "stagnation_inlet_temperature",
    "Stagnation inlet temperature");

  return params;
}

CNSFVStagnationInletBoundaryFlux::CNSFVStagnationInletBoundaryFlux(
    const InputParameters & parameters)
  : BoundaryFluxBase(parameters),
    _bc_uo(getUserObject<BCUserObject>("bc_uo")),
    _fp(getUserObject<SinglePhaseFluidProperties>("fluid_properties")),
    _stag_inlet_pres(getParam<Real>("stagnation_inlet_pressure")),
    _stag_inlet_temp(getParam<Real>("stagnation_inlet_temperature"))
{
}

CNSFVStagnationInletBoundaryFlux::~CNSFVStagnationInletBoundaryFlux() {}

void
CNSFVStagnationInletBoundaryFlux::calcFlux(unsigned int iside,
                                           dof_id_type ielem,
                                           const std::vector<Real> & uvec1,
                                           const RealVectorValue & dwave,
                                           std::vector<Real> & flux) const
{
  std::vector<Real> U2(5, 0.);

  U2 = _bc_uo.getGhostCellValue(iside, ielem, uvec1, dwave);

  Real rho2  = U2[0];
  Real rhou2 = U2[1];
  Real rhov2 = U2[2];
  Real rhow2 = U2[3];
  Real rhoe2 = U2[4];

  Real v2 = 1. / rho2;
  Real uadv2 = v2 * rhou2;
  Real vadv2 = v2 * rhov2;
  Real wadv2 = v2 * rhow2;
  Real vdov2 = uadv2 * uadv2 + vadv2 * vadv2 + wadv2 * wadv2;
  Real e2 = v2 * rhoe2 - 0.5 * vdov2;
  Real pres2 = _fp.pressure(v2, e2);
  Real vdon2 = uadv2 * dwave(0) + vadv2 * dwave(1) + wadv2 * dwave(2);

  flux.resize(5);

  flux[0] = vdon2 * rho2;
  flux[1] = vdon2 * rho2 * uadv2 + pres2 * dwave(0);
  flux[2] = vdon2 * rho2 * vadv2 + pres2 * dwave(1);
  flux[3] = vdon2 * rho2 * wadv2 + pres2 * dwave(2);
  flux[4] = vdon2 * (rhoe2 + pres2);
}

void
CNSFVStagnationInletBoundaryFlux::calcJacobian(unsigned int /*iside*/,
                                               dof_id_type /*ielem*/,
                                               const std::vector<Real> & uvec1,
                                               const RealVectorValue & dwave,
                                               DenseMatrix<Real> & jac1) const
{
  jac1.resize(5, 5);

  const Real rho1  = uvec1[0];
  const Real uadv1 = uvec1[1] / rho1;
  const Real vadv1 = uvec1[2] / rho1;
  const Real wadv1 = uvec1[3] / rho1;
  const Real vdov1 = uadv1 * uadv1 + vadv1 * vadv1 + wadv1 * wadv1;

  const Real nx = dwave(0);
  const Real ny = dwave(1);
  const Real nz = dwave(2);

  const Real p0 = _stag_inlet_pres;
  const Real T0 = _stag_inlet_temp;
  Real rho0, e0;
  _fp.rho_e(p0, T0, rho0, e0);
  const Real v0 = 1.0 / rho0;
  const Real s0 = _fp.s(v0, e0);
  const Real h0 = _fp.h(p0, T0);

  // assume the boundary velocity is exterpolated from the inner cell
  const Real uadvb = uadv1;
  const Real vadvb = vadv1;
  const Real wadvb = wadv1;
  const Real vdovb = uadvb * uadvb + vadvb * vadvb + wadvb * wadvb;
  const Real vdonb = uadvb * nx + vadvb * ny + wadvb * nz;

  const Real sb = s0;
  const Real hb = h0 - 0.5 * vdovb;
  const Real pb = _fp.p_from_h_s(hb, sb);

  Real rhob, drhob_dpb, drhob_dsb, eb, deb_dpb, deb_dsb;
  _fp.rho_e_dps(pb, sb, rhob, drhob_dpb, drhob_dsb, eb, deb_dpb, deb_dsb);
  const Real vb = 1. / rhob;
  const Real Eb = eb + 0.5 * vdovb;
  const Real gamma = _fp.gamma(vb, eb);
  const Real gamm1 = gamma - 1.;
  const Real gamm2 = 2. - gamma;
  const Real rq05b = 0.5 * gamm1 * vdovb;
  const Real enthb = Eb + pb / rhob;

  // ====================
  // Assemble the matrix:
  //
  //   d(Fdotn)
  //   ________
  //
  //    d(U_b)
  // ====================

  DenseMatrix<Real> dFn_dUb(5, 5);

  dFn_dUb(0, 0) = 0.;
  dFn_dUb(0, 1) = nx;
  dFn_dUb(0, 2) = ny;
  dFn_dUb(0, 3) = nz;
  dFn_dUb(0, 4) = 0.;

  dFn_dUb(1, 0) = rq05b * nx - uadvb * vdonb;
  dFn_dUb(1, 1) = gamm2 * nx * uadvb + vdonb;
  dFn_dUb(1, 2) = ny * uadvb - vadvb * gamm1 * nx;
  dFn_dUb(1, 3) = nz * uadvb - wadvb * gamm1 * nx;
  dFn_dUb(1, 4) = gamm1 * nx;

  dFn_dUb(2, 0) = rq05b * ny - vadvb * vdonb;
  dFn_dUb(2, 1) = nx * vadvb - uadvb * gamm1 * ny;
  dFn_dUb(2, 2) = gamm2 * ny * vadvb + vdonb;
  dFn_dUb(2, 3) = nz * vadvb - wadvb * gamm1 * ny;
  dFn_dUb(2, 4) = gamm1 * ny;

  dFn_dUb(3, 0) = rq05b * nz - wadvb * vdonb;
  dFn_dUb(3, 1) = nx * wadvb - uadvb * gamm1 * nz;
  dFn_dUb(3, 2) = ny * wadvb - vadvb * gamm1 * nz;
  dFn_dUb(3, 3) = gamm2 * nz * wadvb + vdonb;
  dFn_dUb(3, 4) = gamm1 * nz;

  dFn_dUb(4, 0) = (rq05b - enthb) * vdonb;
  dFn_dUb(4, 1) = nx * enthb - gamm1 * uadvb * vdonb;
  dFn_dUb(4, 2) = ny * enthb - gamm1 * vadvb * vdonb;
  dFn_dUb(4, 3) = nz * enthb - gamm1 * wadvb * vdonb;
  dFn_dUb(4, 4) = gamma * vdonb;

  // ====================
  // Assemble the matrix:
  //
  //    d(U_b)
  //   ________
  //
  //    d(U_1)
  // ====================

  const Real dpb_dhb = _fp.dpdh_from_h_s(hb, sb);
  const Real dhb_drho1 = vdov1 / rho1;
  const Real dhb_drhou1 = - uadv1 / rho1;
  const Real dhb_drhov1 = - vadv1 / rho1;
  const Real dhb_drhow1 = - wadv1 / rho1;

  const Real dEb_drho1  = deb_dpb * dpb_dhb * dhb_drho1  - vdov1 / rho1;
  const Real dEb_drhou1 = deb_dpb * dpb_dhb * dhb_drhou1 + uadv1 / rho1;
  const Real dEb_drhov1 = deb_dpb * dpb_dhb * dhb_drhov1 + vadv1 / rho1;
  const Real dEb_drhow1 = deb_dpb * dpb_dhb * dhb_drhow1 + wadv1 / rho1;

  DenseMatrix<Real> dUb_dU1(5, 5);

  dUb_dU1(0, 0) = drhob_dpb * dpb_dhb * dhb_drho1;
  dUb_dU1(0, 1) = drhob_dpb * dpb_dhb * dhb_drhou1;
  dUb_dU1(0, 2) = drhob_dpb * dpb_dhb * dhb_drhov1;
  dUb_dU1(0, 3) = drhob_dpb * dpb_dhb * dhb_drhow1;
  dUb_dU1(0, 4) = 0.0;

  dUb_dU1(1, 0) = uadv1 * dUb_dU1(0, 0) - rhob / rho1 * uadv1;
  dUb_dU1(1, 1) = uadv1 * dUb_dU1(0, 1) + rhob / rho1;
  dUb_dU1(1, 2) = uadv1 * dUb_dU1(0, 2);
  dUb_dU1(1, 3) = uadv1 * dUb_dU1(0, 3);
  dUb_dU1(1, 4) = 0.0;

  dUb_dU1(2, 0) = vadv1 * dUb_dU1(0, 0) - rhob / rho1 * vadv1;
  dUb_dU1(2, 1) = vadv1 * dUb_dU1(0, 1);
  dUb_dU1(2, 2) = vadv1 * dUb_dU1(0, 2) + rhob / rho1;
  dUb_dU1(2, 3) = vadv1 * dUb_dU1(0, 3);
  dUb_dU1(2, 4) = 0.0;

  dUb_dU1(3, 0) = wadv1 * dUb_dU1(0, 0) - rhob / rho1 * wadv1;
  dUb_dU1(3, 1) = wadv1 * dUb_dU1(0, 1);
  dUb_dU1(3, 2) = wadv1 * dUb_dU1(0, 2);
  dUb_dU1(3, 3) = wadv1 * dUb_dU1(0, 3) + rhob / rho1;
  dUb_dU1(3, 4) = 0.0;

  dUb_dU1(4, 0) = Eb                 +        dEb_drho1;
  dUb_dU1(4, 1) = Eb * dUb_dU1(0, 1) + rhob * dEb_drhou1;
  dUb_dU1(4, 2) = Eb * dUb_dU1(0, 2) + rhob * dEb_drhov1;
  dUb_dU1(4, 3) = Eb * dUb_dU1(0, 3) + rhob * dEb_drhow1;
  dUb_dU1(4, 4) = 0.0;

  // ================================
  // Follow the chain rule:
  //
  //   d(Fdotn)   d(Fdotn)    d(U_b)
  //   ________ = ________ * ________
  //
  //    d(U_1)     d(U_b)     d(U_1)
  // ================================

  for (unsigned int i = 0; i < 5; i++)
    for (unsigned int j = 0; j < 5; j++)
    {
      jac1(i, j) = 0.;
      for (unsigned int k = 0; k < 5; k++)
        jac1(i, j) += dFn_dUb(i, k) * dUb_dU1(k, j);
    }
}
