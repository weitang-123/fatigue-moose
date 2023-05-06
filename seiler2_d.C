//
//

#include "seiler2_d.h"
#include <iostream>
using namespace std;
registerMooseObject("seiler2App",seiler2_d);
template <>
InputParameters validParams<seiler2_d>()
{
    InputParameters params=validParams<TimeKernel>();
    
    params.addRequiredParam<Real>("G","critical energy release rate");
    params.addRequiredParam<Real>("l0","length");
    params.addRequiredCoupledVar("u1","displacement ux");
    params.addRequiredCoupledVar("u2","displacement uy"); 
    params.addRequiredParam<std::string>("pro_pc","get an initial crack");

    return params;
}

seiler2_d::seiler2_d(const InputParameters &parameters)
        :TimeKernel(parameters),
         _G(getParam<Real>("G")),
         _l0(getParam<Real>("l0")),
         _couple_u1(coupledValue("u1")),
         _couple_u2(coupledValue("u2")),
         _u1_var(coupled("u1")),
         _u2_var(coupled("u2")),
         _pc(getMaterialProperty<Real>(parameters.get<std::string>("pro_pc"))),
         _hist(getMaterialProperty<Real>("hist")),
         _sigma1(getMaterialProperty<Real>("sigma1")),
         _sigma2(getMaterialProperty<Real>("sigma2")),
         _sigma3(getMaterialProperty<Real>("sigma3")),
         _sigmap1(getMaterialProperty<Real>("sigmap1")),
         _sigmap2(getMaterialProperty<Real>("sigmap2")),
         _sigmap3(getMaterialProperty<Real>("sigmap3")),
         _sigmapfc1(getMaterialProperty<Real>("sigmapfc1")),
         _sigmapfc2(getMaterialProperty<Real>("sigmapfc2")),
         _sigmapfc3(getMaterialProperty<Real>("sigmapfc3")),
         _alpha(getMaterialProperty<Real>("alpha"))
{}

Real seiler2_d::computeQpResidual()
{

    // compute R_d
    return 
   (_alpha[_qp]*_u[_qp]-2.0*(1.0-_u[_qp])*_l0/_G*_hist[_qp])*_test[_i][_qp]+_alpha[_qp]*_l0*_l0*_grad_u[_qp]*_grad_test[_i][_qp]- 2.0*_pc[_qp]*(1.0-_u[_qp])*_test[_i][_qp];
}

Real seiler2_d::computeQpJacobian()
{
    // compute K_d_d
    
    return (2.0*_l0/_G*_hist[_qp]+_alpha[_qp]+2.0*_pc[_qp])*_test[_i][_qp]*_phi[_j][_qp]+_alpha[_qp]*_l0*_l0*_grad_phi[_j][_qp]*_grad_test[_i][_qp] ;
}

Real seiler2_d::computeQpOffDiagJacobian(unsigned int jvar)
{
    if(jvar==_u1_var)
    {
         // K_d_u1
         return -2.0*(1.0-_u[_qp])*_l0/_G*(_sigmapfc1[_qp]*_grad_phi[_j][_qp](0)+_sigmapfc3[_qp]*_grad_phi[_j][_qp](1))*_test[_i][_qp] ;
    }

    if(jvar==_u2_var)
    {
         // K_d_u2
        return -2.0*(1.0-_u[_qp])*_l0/_G*(_sigmapfc2[_qp]*_grad_phi[_j][_qp](1)+_sigmapfc3[_qp]*_grad_phi[_j][_qp](0))*_test[_i][_qp] ;
    }
    
    return 0.0; // if not the Kcmu ,return zero contribution
}
    

