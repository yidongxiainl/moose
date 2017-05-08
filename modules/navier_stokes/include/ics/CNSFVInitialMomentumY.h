/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CNSFVINITIALMOMENTUMY_H
#define CNSFVINITIALMOMENTUMY_H

#include "CNSFVICBase.h"

class CNSFVInitialMomentumY;

template <>
InputParameters validParams<CNSFVInitialMomentumY>();

/**
 * CNSFVInitialMomentumY sets intial constant values for y-momentum given the:
 * .) Initial pressure
 * .) Initial temperature
 * .) Initial velocity
 * and a FluidProperties UserObject.
 */
class CNSFVInitialMomentumY : public CNSFVICBase
{
public:
  CNSFVInitialMomentumY(const InputParameters & parameters);

  virtual Real value(const Point & p);

protected:
};

#endif
