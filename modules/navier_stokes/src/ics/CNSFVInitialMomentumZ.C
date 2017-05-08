/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CNSFVInitialMomentumZ.h"

template <>
InputParameters
validParams<CNSFVInitialMomentumZ>()
{
  InputParameters params = validParams<CNSFVICBase>();

  params.addClassDescription(
    "CNSFVInitialMomentumZ sets intial constant values for z-momentum.");

  return params;
}

CNSFVInitialMomentumZ::CNSFVInitialMomentumZ(const InputParameters & parameters)
  : CNSFVICBase(parameters)
{
}

Real
CNSFVInitialMomentumZ::value(const Point & /*p*/)
{
  return _fp.rho(_init_pres, _init_temp) * _init_velo(2);
}
