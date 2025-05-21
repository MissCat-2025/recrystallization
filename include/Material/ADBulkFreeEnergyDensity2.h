// ADBulkFreeEnergyDensity2.h
#pragma once

#include "ADMaterial.h"

/**
 * 使用自动微分计算多晶与气泡相的体自由能密度材料
 */
class ADBulkFreeEnergyDensity2 : public ADMaterial
{
public:
  static InputParameters validParams();
  ADBulkFreeEnergyDensity2(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  
  // 耦合的气体浓度变量
  const ADVariableValue & _c_g;
  
  // 晶粒序参数
  const unsigned int _num_grains;
  std::vector<const ADVariableValue *> _grain_vars;
  
  // 气泡序参数
  const unsigned int _num_bubbles;
  std::vector<const ADVariableValue *> _bubble_vars;
  
  // 材料参数
  const Real _A;
  const Real _C_p;
  const Real _C_q;
  const Real _S_I;
  const Real _c_m_e;
  const Real _c_b_e;
  
  // 材料属性
  ADMaterialProperty<Real> & _f_bulk;
};