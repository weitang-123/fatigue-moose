#include "AuxKernel.h"

// Forward Declarations
class seiler2_aux;

template <>
InputParameters validParams<seiler2_aux>();

/**
 * Coupled auxiliary value
 */
class seiler2_aux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  seiler2_aux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const unsigned int _i;
  const MaterialProperty<Real> &_epsilon1,&_epsilon2,&_epsilon3,&_sigma1,&_sigma2,&_sigma3,&_Pswt1,&_Pswt2,&_Pswt3,&_D,&_sigmabar,&_sigmamax,&_sigmareva,&_alpha,&_hist,&_Ni;

};
