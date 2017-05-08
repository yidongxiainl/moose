/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CNSFVInitialTotalEnergy.h"

template <>
InputParameters
validParams<CNSFVInitialTotalEnergy>()
{
  InputParameters params = validParams<CNSFVICBase>();

  params.addClassDescription(
    "CNSFVInitialTotalEnergy sets intial constant values for total energy.");

  return params;
}

CNSFVInitialTotalEnergy::CNSFVInitialTotalEnergy(const InputParameters & parameters)
  : CNSFVICBase(parameters)
{
}

Real
CNSFVInitialTotalEnergy::value(const Point & /*p*/)
{
  const Real rho = _fp.rho(_init_pres, _init_temp);
  return rho * (_fp.e(_init_pres, rho) + 0.5 * _init_velo.norm_sq());
}
