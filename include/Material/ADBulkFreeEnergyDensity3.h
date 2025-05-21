// ADBulkFreeEnergyDensity3.h
#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"

/**
 * 材料类用于计算多晶与气泡相的体自由能密度
 * 使用DerivativeMaterialInterface实现自动导数计算
 * 输入参数包括晶粒参数η_i(i=1,2,...,p)和气泡参数η_j(j=p+1,p+2,...,q)
 * 以及气体浓度c_g
 */
class ADBulkFreeEnergyDensity3 : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ADBulkFreeEnergyDensity3(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// 气体浓度变量
  const VariableValue & _c_g;
  const VariableName _c_g_name;

  /// 晶粒序参数（多晶结构相）
  unsigned int _num_grains;
  std::vector<const VariableValue *> _grain_vars;
  std::vector<VariableName> _grain_var_names;

  /// 气泡序参数
  unsigned int _num_bubbles;
  std::vector<const VariableValue *> _bubble_vars;
  std::vector<VariableName> _bubble_var_names;

  /// 体自由能相关系数
  const Real _A;
  const Real _C_p;
  const Real _C_q;
  const Real _S_I;
  const Real _c_m_e;
  const Real _c_b_e;

  /// 材料属性：体自由能密度
  MaterialProperty<Real> & _f_bulk;
};