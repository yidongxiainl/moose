/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CNSFVINITIALDENSITY_H
#define CNSFVINITIALDENSITY_H

#include "CNSFVICBase.h"

class CNSFVInitialDensity;

template <>
InputParameters validParams<CNSFVInitialDensity>();

/**
 * CNSFVInitialDensity sets intial constant values for density given the:
 * .) Initial pressure
 * .) Initial temperature
 * .) Initial velocity
 * and a FluidProperties UserObject.
 */
class CNSFVInitialDensity : public CNSFVICBase
{
public:
  CNSFVInitialDensity(const InputParameters & parameters);

  virtual Real value(const Point & p);

protected:
};

#endif
