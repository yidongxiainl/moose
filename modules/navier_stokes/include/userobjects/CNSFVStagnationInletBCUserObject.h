/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef CNSFVSTAGNATIONINLECTBCUSEROBJECT_H
#define CNSFVSTAGNATIONINLECTBCUSEROBJECT_H

#include "BCUserObject.h"
#include "SinglePhaseFluidProperties.h"

class CNSFVStagnationInletBCUserObject;

template <>
InputParameters validParams<CNSFVStagnationInletBCUserObject>();

/**
 * A user object that computes the ghost cell values based on
 * the stagnation inlet boundary condition
 *
 * Reference:
 *
 * Berry, Ray A., Richard Saurel, and Olivier LeMetayer.
 * "The discrete equation method (DEM) for fully compressible,
 * two-phase flows in ducts of spatially varying cross-section."
 * Nuclear Engineering and Design 240.11 (2010): 3797-3818.
 */
class CNSFVStagnationInletBCUserObject : public BCUserObject
{
public:
  CNSFVStagnationInletBCUserObject(const InputParameters & parameters);

  virtual std::vector<Real> getGhostCellValue(unsigned int iside,
                                              dof_id_type ielem,
                                              const std::vector<Real> & uvec1,
                                              const RealVectorValue & dwave) const;

protected:
  const SinglePhaseFluidProperties & _fp;

  Real _stag_inlet_pres;
  Real _stag_inlet_temp;
};

#endif
