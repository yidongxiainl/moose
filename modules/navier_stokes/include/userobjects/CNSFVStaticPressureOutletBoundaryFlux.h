/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef CNSFVSTATICPRESSUREOUTLETBOUNDARYFLUX_H
#define CNSFVSTATICPRESSUREOUTLETBOUNDARYFLUX_H

#include "BoundaryFluxBase.h"
#include "BCUserObject.h"
#include "SinglePhaseFluidProperties.h"

// Forward Declarations
class CNSFVStaticPressureOutletBoundaryFlux;

template <>
InputParameters validParams<CNSFVStaticPressureOutletBoundaryFlux>();

/**
 * A user objec that computes the Riemann-invariant boundary flux
 */
class CNSFVStaticPressureOutletBoundaryFlux : public BoundaryFluxBase
{
public:
  CNSFVStaticPressureOutletBoundaryFlux(const InputParameters & parameters);
  virtual ~CNSFVStaticPressureOutletBoundaryFlux();

  virtual void calcFlux(unsigned int iside,
                        dof_id_type ielem,
                        const std::vector<Real> & uvec1,
                        const RealVectorValue & dwave,
                        std::vector<Real> & flux) const;

  virtual void calcJacobian(unsigned int iside,
                            dof_id_type ielem,
                            const std::vector<Real> & uvec1,
                            const RealVectorValue & dwave,
                            DenseMatrix<Real> & jac1) const;

protected:
  const BCUserObject & _bc_uo;
  const SinglePhaseFluidProperties & _fp;

  Real _static_out_pres;
};

#endif
