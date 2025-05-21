// ADDislocEnergyDensity2.h
#pragma once

#include "ADMaterial.h"

/**
 * 使用自动微分计算位错能量密度
 */
class ADDislocEnergyDensity2 : public ADMaterial
{
public:
  static InputParameters validParams();
  ADDislocEnergyDensity2(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  
  // 计算位错密度
  Real computeDislocationDensity() const;
  
  // 材料参数
  const Real _shear_modulus;
  const Real _burgers_vector;
  const Real _rho_d_0;
  const Real _lambda;
  const Real _fission_rate;
  
  // 相场变量
  const unsigned int _num_etas;
  std::vector<const ADVariableValue *> _eta_values;
  
  // 材料属性
  MaterialProperty<Real> & _dislocation_density;
  ADMaterialProperty<Real> & _f_stored;
  
  // 时间
  const Real & _t;
};