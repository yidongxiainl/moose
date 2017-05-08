/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CNSFVINITIALTOTALENERGY_H
#define CNSFVINITIALTOTALENERGY_H

#include "CNSFVICBase.h"

class CNSFVInitialTotalEnergy;

template <>
InputParameters validParams<CNSFVInitialTotalEnergy>();

/**
 * CNSFVInitialTotalEnergy sets intial constant values for total energy given the:
 * .) Initial pressure
 * .) Initial temperature
 * .) Initial velocity
 * and a FluidProperties UserObject.
 */
class CNSFVInitialTotalEnergy : public CNSFVICBase
{
public:
  CNSFVInitialTotalEnergy(const InputParameters & parameters);

  virtual Real value(const Point & p);

protected:
};

#endif
