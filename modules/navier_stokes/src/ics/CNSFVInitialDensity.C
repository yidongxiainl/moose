/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CNSFVInitialDensity.h"

template <>
InputParameters
validParams<CNSFVInitialDensity>()
{
  InputParameters params = validParams<CNSFVICBase>();

  params.addClassDescription(
    "CNSFVInitialDensity sets intial constant values for density.");

  return params;
}

CNSFVInitialDensity::CNSFVInitialDensity(const InputParameters & parameters)
  : CNSFVICBase(parameters)
{
}

Real
CNSFVInitialDensity::value(const Point & /*p*/)
{
  return _fp.rho(_init_pres, _init_temp);
}
