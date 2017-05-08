/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CNSFVICBASE_H
#define CNSFVICBASE_H

#include "InitialCondition.h"
#include "SinglePhaseFluidProperties.h"

class CNSFVICBase;

template <>
InputParameters validParams<CNSFVICBase>();

/**
 * CNSFVICBase sets intial constant values for all variables given the:
 * .) Initial pressure
 * .) Initial temperature
 * .) Initial velocity
 * and a FluidProperties UserObject.
 */
class CNSFVICBase : public InitialCondition
{
public:
  CNSFVICBase(const InputParameters & parameters);

  virtual Real value(const Point & p) = 0;

protected:
  Real _init_pres;
  Real _init_temp;
  RealVectorValue _init_velo;

  const SinglePhaseFluidProperties & _fp;
};

#endif
