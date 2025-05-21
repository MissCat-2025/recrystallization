// ADDislocEnergyDensity2.C
#include "ADDislocEnergyDensity2.h"

registerMooseObject("RecrystallizationApp", ADDislocEnergyDensity2);

InputParameters
ADDislocEnergyDensity2::validParams()
{
  InputParameters params = ADMaterial::validParams();
  
  params.addParam<Real>("shear_modulus", 36.0e9, "剪切模量 G (Pa)");
  params.addParam<Real>("burgers_vector", 3.42e-10, "伯格斯矢量大小 b_v (m)");
  params.addParam<Real>("rho_d_0", 1.0e13, "初始位错密度 ρ_d(0) (m^-2)");
  params.addParam<Real>("lambda", 1.0e-20, "位错累积速率 λ");
  params.addParam<Real>("fission_rate", 5.0e20, "裂变速率 ḟ (m^-3 s^-1)");
  
  params.addRequiredCoupledVar("etas", "晶体变量η的列表");
  
  params.addParam<MaterialPropertyName>("dislocation_density_name", "dislocation_density", "位错密度材料属性名称");
  params.addParam<MaterialPropertyName>("f_name", "f_stored", "储存弹性能量密度材料属性名称");
  
  params.addClassDescription("使用自动微分计算位错能量密度");
  
  return params;
}

ADDislocEnergyDensity2::ADDislocEnergyDensity2(const InputParameters & parameters) :
  ADMaterial(parameters),
  _shear_modulus(getParam<Real>("shear_modulus")),
  _burgers_vector(getParam<Real>("burgers_vector")),
  _rho_d_0(getParam<Real>("rho_d_0")),
  _lambda(getParam<Real>("lambda")),
  _fission_rate(getParam<Real>("fission_rate")),
  _num_etas(coupledComponents("etas")),
  _eta_values(_num_etas),
  _dislocation_density(declareProperty<Real>(getParam<MaterialPropertyName>("dislocation_density_name"))),
  _f_stored(declareADProperty<Real>(getParam<MaterialPropertyName>("f_name"))),
  _t(_fe_problem.time())
{
  // 获取所有相场变量
  for (unsigned int i = 0; i < _num_etas; ++i)
    _eta_values[i] = &adCoupledValue("etas", i);
}

void
ADDislocEnergyDensity2::computeQpProperties()
{
  // 计算位错密度
  Real rho_d = computeDislocationDensity();
  _dislocation_density[_qp] = rho_d;
  
  // 计算常数因子
  Real factor = 0.5 * _shear_modulus * _burgers_vector * _burgers_vector * rho_d;
  
  // 计算总储存能量密度
  _f_stored[_qp] = 0.0;
  for (unsigned int i = 0; i < _num_etas; ++i)
  {
    ADReal eta_i = (*_eta_values[i])[_qp];
    _f_stored[_qp] += factor * eta_i * eta_i;
  }
}

Real
ADDislocEnergyDensity2::computeDislocationDensity() const
{
  // 根据公式计算位错密度：ρ_d(t) = ρ_d(0)⋅exp[λḟ⋅t]，后面可拓展为与温度相关
  return _rho_d_0 * std::exp(_lambda * _fission_rate * _t);
}