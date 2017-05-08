/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CNSFVInitialMomentumY.h"

template <>
InputParameters
validParams<CNSFVInitialMomentumY>()
{
  InputParameters params = validParams<CNSFVICBase>();

  params.addClassDescription(
    "CNSFVInitialMomentumY sets intial constant values for y-momentum.");

  return params;
}

CNSFVInitialMomentumY::CNSFVInitialMomentumY(const InputParameters & parameters)
  : CNSFVICBase(parameters)
{
}

Real
CNSFVInitialMomentumY::value(const Point & /*p*/)
{
  return _fp.rho(_init_pres, _init_temp) * _init_velo(1);
}
