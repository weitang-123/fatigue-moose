#include "seiler2_aux.h"

registerMooseObject("seiler2App", seiler2_aux);

template <>
InputParameters
validParams<seiler2_aux>()
{
    InputParameters params = validParams<AuxKernel>();
    params.addRequiredParam<unsigned int>("i", "The i: 1-12 indicates seiler2 (1-6) and strain (7-12)");
    
    return params;
}

seiler2_aux::seiler2_aux(const InputParameters & parameters)
  : AuxKernel(parameters),
      _i(getParam<unsigned int>("i")),
      _epsilon1(getMaterialProperty<Real>("epsilon1")),
      _epsilon2(getMaterialProperty<Real>("epsilon2")),
      _epsilon3(getMaterialProperty<Real>("epsilon3")),
      _sigma1(getMaterialProperty<Real>("sigma1")),
      _sigma2(getMaterialProperty<Real>("sigma2")),
      _sigma3(getMaterialProperty<Real>("sigma3")),
      _Pswt1(getMaterialProperty<Real>("Pswt1")),
      _Pswt2(getMaterialProperty<Real>("Pswt2")),
      _Pswt3(getMaterialProperty<Real>("Pswt3")),
      _D(getMaterialProperty<Real>("D")),
      _sigmabar(getMaterialProperty<Real>("sigmabar")),
      _sigmamax(getMaterialProperty<Real>("sigmamax")),
      _sigmareva(getMaterialProperty<Real>("sigmareva")),
      _alpha(getMaterialProperty<Real>("alpha")),
      _hist(getMaterialProperty<Real>("hist")),
      _Ni(getMaterialProperty<Real>("Ni"))
      
{
}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */
Real
seiler2_aux::computeValue()
{
    float val[20]; // array must use [number of elements]
    
    val[0]=_epsilon1[_qp]; val[1]=_epsilon2[_qp]; val[2]=_epsilon3[_qp];   
   
    val[6]=_sigma1[_qp]; val[7]=_sigma2[_qp]; val[8]=_sigma3[_qp];
    
    val[9]=_D[_qp]; val[10]=_sigmabar[_qp]; val[11]=_alpha[_qp];
    
    val[12]=_sigmamax[_qp];val[13]=_Ni[_qp];val[14]=_sigmareva[_qp];val[15]=_hist[_qp]; 
    
    val[16]=_Pswt1[_qp]; val[17]=_Pswt2[_qp]; val[18]=_Pswt3[_qp];
    
    return val[_i];
}


