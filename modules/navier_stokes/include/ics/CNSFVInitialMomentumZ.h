/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CNSFVINITIALMOMENTUMZ_H
#define CNSFVINITIALMOMENTUMZ_H

#include "CNSFVICBase.h"

class CNSFVInitialMomentumZ;

template <>
InputParameters validParams<CNSFVInitialMomentumZ>();

/**
 * CNSFVInitialMomentumZ sets intial constant values for z-momentum given the:
 * .) Initial pressure
 * .) Initial temperature
 * .) Initial velocity
 * and a FluidProperties UserObject.
 */
class CNSFVInitialMomentumZ : public CNSFVICBase
{
public:
  CNSFVInitialMomentumZ(const InputParameters & parameters);

  virtual Real value(const Point & p);

protected:
};

#endif
