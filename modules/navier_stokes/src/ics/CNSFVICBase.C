/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CNSFVICBase.h"

template <>
InputParameters
validParams<CNSFVICBase>()
{
  InputParameters params = validParams<InitialCondition>();

  params.addClassDescription(
    "CNSFVICBase sets intial constant values for conserved variables.");

  params.addRequiredParam<Real>(
    "initial_pressure",
    "The initial pressure, assumed constant everywhere");

  params.addRequiredParam<Real>(
    "initial_temperature",
    "The initial temperature, assumed constant everywhere");

  params.addRequiredParam<RealVectorValue>(
    "initial_velocity",
    "The initial velocity, assumed constant everywhere");

  params.addRequiredParam<UserObjectName>(
    "fluid_properties",
    "The name of the user object for fluid properties");

  return params;
}

CNSFVICBase::CNSFVICBase(const InputParameters & parameters)
  : InitialCondition(parameters),
    _init_pres(getParam<Real>("initial_pressure")),
    _init_temp(getParam<Real>("initial_temperature")),
    _init_velo(getParam<RealVectorValue>("initial_velocity")),
    _fp(getUserObject<SinglePhaseFluidProperties>("fluid_properties"))
{
}
