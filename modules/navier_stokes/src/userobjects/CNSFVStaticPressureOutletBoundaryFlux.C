/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CNSFVStaticPressureOutletBoundaryFlux.h"

template <>
InputParameters
validParams<CNSFVStaticPressureOutletBoundaryFlux>()
{
  InputParameters params = validParams<BoundaryFluxBase>();

  params.addClassDescription(
    "A user objec that computes the static pressure outlet boundary flux.");

  params.addRequiredParam<UserObjectName>(
    "bc_uo",
    "Name for boundary condition user object");

  params.addRequiredParam<UserObjectName>(
    "fluid_properties",
    "Name for fluid properties user object");

  params.addRequiredParam<Real>(
    "static_outlet_pressure",
    "Static outlet pressure");

  return params;
}

CNSFVStaticPressureOutletBoundaryFlux::CNSFVStaticPressureOutletBoundaryFlux(
    const InputParameters & parameters)
  : BoundaryFluxBase(parameters),
    _bc_uo(getUserObject<BCUserObject>("bc_uo")),
    _fp(getUserObject<SinglePhaseFluidProperties>("fluid_properties")),
    _static_out_pres(getParam<Real>("static_outlet_pressure"))
{
}

CNSFVStaticPressureOutletBoundaryFlux::~CNSFVStaticPressureOutletBoundaryFlux() {}

void
CNSFVStaticPressureOutletBoundaryFlux::calcFlux(unsigned int iside,
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
CNSFVStaticPressureOutletBoundaryFlux::calcJacobian(unsigned int /*iside*/,
                                                    dof_id_type /*ielem*/,
                                                    const std::vector<Real> & uvec1,
                                                    const RealVectorValue & dwave,
                                                    DenseMatrix<Real> & jac1) const
{
  jac1.resize(5, 5);

  const Real rho1  = uvec1[0];
  const Real uadv1 = uvec1[1] / uvec1[0];
  const Real vadv1 = uvec1[2] / uvec1[0];
  const Real wadv1 = uvec1[3] / uvec1[0];
  const Real rhoe1 = uvec1[4];

  const Real nx = dwave(0);
  const Real ny = dwave(1);
  const Real nz = dwave(2);

  const Real vdon1 = uadv1 * nx + vadv1 * ny + wadv1 * nz;
  const Real enth1 = (rhoe1 + _static_out_pres) / rho1;

  jac1(0, 0) = 0.;
  jac1(0, 1) = nx;
  jac1(0, 2) = ny;
  jac1(0, 3) = nz;
  jac1(0, 4) = 0.;

  jac1(1, 0) = -vdon1 * uadv1;
  jac1(1, 1) = nx * uadv1 + vdon1;
  jac1(1, 2) = ny * uadv1;
  jac1(1, 3) = nz * uadv1;
  jac1(1, 4) = 0.;

  jac1(2, 0) = -vdon1 * vadv1;
  jac1(2, 1) = nx * vadv1;
  jac1(2, 2) = ny * vadv1 + vdon1;
  jac1(2, 3) = nz * vadv1;
  jac1(2, 4) = 0.;

  jac1(3, 0) = -vdon1 * wadv1;
  jac1(3, 1) = nx * wadv1;
  jac1(3, 2) = ny * wadv1;
  jac1(3, 3) = nz * wadv1 + vdon1;
  jac1(3, 4) = 0.;

  jac1(4, 0) = -vdon1 * enth1;
  jac1(4, 1) = nx * enth1;
  jac1(4, 2) = ny * enth1;
  jac1(4, 3) = nz * enth1;
  jac1(4, 4) = vdon1;
}
