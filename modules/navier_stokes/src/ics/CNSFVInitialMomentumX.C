/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CNSFVInitialMomentumX.h"

template <>
InputParameters
validParams<CNSFVInitialMomentumX>()
{
  InputParameters params = validParams<CNSFVICBase>();

  params.addClassDescription(
    "CNSFVInitialMomentumX sets intial constant values for x-momentum.");

  return params;
}

CNSFVInitialMomentumX::CNSFVInitialMomentumX(const InputParameters & parameters)
  : CNSFVICBase(parameters)
{
}

Real
CNSFVInitialMomentumX::value(const Point & /*p*/)
{
  return _fp.rho(_init_pres, _init_temp) * _init_velo(0);
}
