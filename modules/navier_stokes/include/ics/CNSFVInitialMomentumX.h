/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CNSFVINITIALMOMENTUMX_H
#define CNSFVINITIALMOMENTUMX_H

#include "CNSFVICBase.h"

class CNSFVInitialMomentumX;

template <>
InputParameters validParams<CNSFVInitialMomentumX>();

/**
 * CNSFVInitialMomentumX sets intial constant values for x-momentum given the:
 * .) Initial pressure
 * .) Initial temperature
 * .) Initial velocity
 * and a FluidProperties UserObject.
 */
class CNSFVInitialMomentumX : public CNSFVICBase
{
public:
  CNSFVInitialMomentumX(const InputParameters & parameters);

  virtual Real value(const Point & p);

protected:
};

#endif
