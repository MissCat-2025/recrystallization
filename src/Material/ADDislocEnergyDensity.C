// ADDislocEnergyDensity.C
#include "ADDislocEnergyDensity.h"

registerMooseObject("RecrystallizationApp", ADDislocEnergyDensity);

InputParameters
ADDislocEnergyDensity::validParams()
{
  InputParameters params = Material::validParams();
  
  params.addParam<Real>("shear_modulus", 36.0e9, "剪切模量 G (Pa)");
  params.addParam<Real>("burgers_vector", 3.42e-10, "伯格斯矢量大小 b_v (m)");
  params.addParam<Real>("rho_d_0", 1.0e13, "初始位错密度 ρ_d(0) (m^-2)");
  params.addParam<Real>("lambda", 1.0e-20, "位错累积速率 λ");
  params.addParam<Real>("fission_rate", 5.0e20, "裂变速率 ḟ (m^-3 s^-1)");
  
  params.addRequiredCoupledVar("etas", "相场变量η的列表");
  
  params.addParam<MaterialPropertyName>("dislocation_density_name", "dislocation_density", "位错密度材料属性名称");
  params.addParam<MaterialPropertyName>("f_name", "f_stored", "储存弹性能量密度材料属性名称");
  
  return params;
}

ADDislocEnergyDensity::ADDislocEnergyDensity(const InputParameters & parameters) :
  DerivativeMaterialInterface<Material>(parameters),
  _shear_modulus(getParam<Real>("shear_modulus")),
  _burgers_vector(getParam<Real>("burgers_vector")),
  _rho_d_0(getParam<Real>("rho_d_0")),
  _lambda(getParam<Real>("lambda")),
  _fission_rate(getParam<Real>("fission_rate")),
  _num_etas(coupledComponents("etas")),
  _eta_values(_num_etas),
  _dislocation_density(declareProperty<Real>(getParam<MaterialPropertyName>("dislocation_density_name"))),
  _f_stored(declareProperty<Real>(getParam<MaterialPropertyName>("f_name"))),
  _df_stored_deta(_num_etas),
  _d2f_stored_deta(_num_etas),
  _t(_fe_problem.time())
{
  // 获取所有相场变量
  for (unsigned int i = 0; i < _num_etas; ++i)
  {
    _eta_values[i] = &coupledValue("etas", i);
    
    // 声明对每个相场变量的一阶导数
    _df_stored_deta[i] = &declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("f_name"), getVar("etas", i)->name());
    
    // 声明对每个相场变量的二阶导数
    _d2f_stored_deta[i].resize(_num_etas);
    for (unsigned int j = 0; j < _num_etas; ++j)
    {
      _d2f_stored_deta[i][j] = &declarePropertyDerivative<Real>(
          getParam<MaterialPropertyName>("f_name"), 
          getVar("etas", i)->name(), 
          getVar("etas", j)->name());
    }
  }
}

void
ADDislocEnergyDensity::computeQpProperties()
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
    Real eta_i = (*_eta_values[i])[_qp];
    _f_stored[_qp] += factor * eta_i * eta_i;
  }
  
  // 计算一阶导数: df/deta_i = G * b_v^2 * rho_d * eta_i
  for (unsigned int i = 0; i < _num_etas; ++i)
  {
    Real eta_i = (*_eta_values[i])[_qp];
    (*_df_stored_deta[i])[_qp] = factor * 2.0 * eta_i;
  }
  
  // 计算二阶导数: d^2f/deta_i deta_j
  for (unsigned int i = 0; i < _num_etas; ++i)
  {
    for (unsigned int j = 0; j < _num_etas; ++j)
    {
      if (i == j)
        // d^2f/deta_i^2 = G * b_v^2 * rho_d * 2
        (*_d2f_stored_deta[i][j])[_qp] = factor * 2.0;
      else
        // d^2f/deta_i deta_j = 0 (i != j)
        (*_d2f_stored_deta[i][j])[_qp] = 0.0;
    }
  }
}

Real
ADDislocEnergyDensity::computeDislocationDensity() const
{
  // 根据公式计算位错密度：ρ_d(t) = ρ_d(0)⋅exp[λḟ⋅t]，后面可拓展为与温度相关
  return _rho_d_0 * std::exp(_lambda * _fission_rate * _t);
}