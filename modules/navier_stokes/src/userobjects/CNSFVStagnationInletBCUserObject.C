/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CNSFVStagnationInletBCUserObject.h"

template <>
InputParameters
validParams<CNSFVStagnationInletBCUserObject>()
{
  InputParameters params = validParams<BCUserObject>();

  params.addClassDescription(
    "A user object that computes conserved variable values in ghost cell "
    "based on the stagnation inlet boundary condition.");

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

CNSFVStagnationInletBCUserObject::CNSFVStagnationInletBCUserObject(
    const InputParameters & parameters)
  : BCUserObject(parameters),
    _fp(getUserObject<SinglePhaseFluidProperties>("fluid_properties")),
    _stag_inlet_pres(getParam<Real>("stagnation_inlet_pressure")),
    _stag_inlet_temp(getParam<Real>("stagnation_inlet_temperature"))
{
}

std::vector<Real>
CNSFVStagnationInletBCUserObject::getGhostCellValue(unsigned int /*iside*/,
                                                    dof_id_type /*ielem*/,
                                                    const std::vector<Real> & uvec1,
                                                    const RealVectorValue & /*dwave*/) const
{
  Real rho1  = uvec1[0];
  Real uadv1 = uvec1[1] / rho1;
  Real vadv1 = uvec1[2] / rho1;
  Real wadv1 = uvec1[3] / rho1;

  Real p0b = _stag_inlet_pres;
  Real t0b = _stag_inlet_temp;
  Real r0b = _fp.rho(p0b, t0b);
  Real e0b = _fp.e(p0b, r0b);
  Real h0b = _fp.h(p0b, t0b);
  Real v0b = 1.0 / r0b;
  Real s0b = _fp.s(v0b, e0b);

  Real sb = s0b;
  Real ub = std::sqrt(uadv1 * uadv1 + vadv1 * vadv1 + wadv1 * wadv1);
  Real hb = h0b - 0.5 * ub * ub;
  Real pb = _fp.p_from_h_s(hb, sb);
  Real rb = 0.0;
  Real eb = 0.0;
  _fp.rho_e_ps(pb, sb, rb, eb);
  Real Eb  = eb + 0.5 * ub * ub;

  std::vector<Real> urigh(5, 0.);

  urigh[0] = rb;
  urigh[1] = rb * uadv1;
  urigh[2] = rb * vadv1;
  urigh[3] = rb * wadv1;
  urigh[4] = rb * Eb;

  return urigh;
}
