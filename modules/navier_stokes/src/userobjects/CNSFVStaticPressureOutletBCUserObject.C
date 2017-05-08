/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CNSFVStaticPressureOutletBCUserObject.h"

template <>
InputParameters
validParams<CNSFVStaticPressureOutletBCUserObject>()
{
  InputParameters params = validParams<BCUserObject>();

  params.addClassDescription(
    "A user object that computes conserved variable values in ghost cell "
    "based on the static pressure outlet boundary condition.");

  params.addRequiredParam<UserObjectName>(
    "fluid_properties",
    "Name for fluid properties user object");

  params.addRequiredParam<Real>(
    "static_outlet_pressure",
    "Static outlet pressure");

  return params;
}

CNSFVStaticPressureOutletBCUserObject::CNSFVStaticPressureOutletBCUserObject(
    const InputParameters & parameters)
  : BCUserObject(parameters),
    _fp(getUserObject<SinglePhaseFluidProperties>("fluid_properties")),
    _static_out_pres(getParam<Real>("static_outlet_pressure"))
{
}

std::vector<Real>
CNSFVStaticPressureOutletBCUserObject::getGhostCellValue(unsigned int /*iside*/,
                                                    dof_id_type /*ielem*/,
                                                    const std::vector<Real> & uvec1,
                                                    const RealVectorValue & /*dwave*/) const
{
  Real rho1 = uvec1[0];
  Real uadv1 = uvec1[1] / rho1;
  Real vadv1 = uvec1[2] / rho1;
  Real wadv1 = uvec1[3] / rho1;

  Real rhoe2 = rho1 * (_fp.e(_static_out_pres, rho1) + 0.5 * (uadv1 * uadv1 + vadv1 * vadv1 + wadv1 * wadv1));

  std::vector<Real> urigh(5, 0.);

  urigh[0] = uvec1[0];
  urigh[1] = uvec1[1];
  urigh[2] = uvec1[2];
  urigh[3] = uvec1[3];
  urigh[4] = rhoe2;

  return urigh;
}
