// ADGasSourceSink.C
#include "ADGasSourceSink.h"
#include "MooseRandom.h"
#include "FEProblem.h"
#include "MooseMesh.h"

registerMooseObject("RecrystallizationApp", ADGasSourceSink);

InputParameters
ADGasSourceSink::validParams()
{
  InputParameters params = ADKernel::validParams();
  
  params.addRequiredParam<Real>("fission_rate", "裂变率 f ̇ (m^(-3) s^(-1))");
  params.addRequiredParam<Real>("atomic_volume", "原子体积 Ω (m^3)");
  params.addRequiredParam<Real>("resolution_rate", "气体分辨率系数 Λ (s^(-1))");
  params.addRequiredCoupledVar("bubble_etas", "气泡相序参数列表");
  params.addParam<unsigned int>("seed", 0, "随机数生成器种子");
  params.addParam<bool>("update_seed_per_timestep", true, "是否在每个时间步更新随机数");
  
  params.addClassDescription("实现气体浓度方程中的源项: G ̇ - R ̇");
  
  return params;
}

ADGasSourceSink::ADGasSourceSink(const InputParameters & parameters)
  : ADKernel(parameters),
    _fission_rate(getParam<Real>("fission_rate")),
    _atomic_volume(getParam<Real>("atomic_volume")),
    _resolution_rate(getParam<Real>("resolution_rate")),
    _num_bubbles(coupledComponents("bubble_etas")),
    _seed(getParam<unsigned int>("seed")),
    _update_seed_per_timestep(getParam<bool>("update_seed_per_timestep"))
{
  // 存储气泡序参数
  _bubble_vars.resize(_num_bubbles);
  for (unsigned int i = 0; i < _num_bubbles; ++i)
    _bubble_vars[i] = &adCoupledValue("bubble_etas", i);
  
  // 修复1: 使用正确的方法名nElem()
  _random_values.resize(_mesh.nElem(), 0.0);
  
  // 设置随机数生成器并生成初始随机值
  MooseRandom::seed(_seed);
  for (unsigned int i = 0; i < _random_values.size(); ++i)
    _random_values[i] = MooseRandom::rand();
}

ADReal
ADGasSourceSink::computeQpResidual()
{
  // 计算G ̇项 (裂变气体产生率)
  // G ̇ = ϖ·Ran, 其中 ϖ = f ̇·Ω
  const Real omega = _fission_rate * _atomic_volume;
  
  // 修复2: 使用_fe_problem.getCurrentTime()和_fe_problem.getTime()
  // 或者使用_t和_t_old来判断时间步是否改变
  if (_update_seed_per_timestep && std::abs(_t - _t_old) > 1e-10)
  {
    for (unsigned int i = 0; i < _random_values.size(); ++i)
      _random_values[i] = MooseRandom::rand();
  }
  
  // 获取当前元素的随机数
  const Real random_value = _random_values[_current_elem->id()];
  
  // 计算气体产生率 G ̇
  const Real gas_generation_rate = omega * random_value;
  
  // 计算R ̇项 (气体分辨率)
  // R ̇ = η_i^2·Λ·c_g, 其中i从p+1到q (气泡相参数)
  ADReal bubble_sum_squared = 0.0;
  for (unsigned int i = 0; i < _num_bubbles; ++i)
    bubble_sum_squared += (*_bubble_vars[i])[_qp] * (*_bubble_vars[i])[_qp];
  
  const ADReal resolution_rate = bubble_sum_squared * _resolution_rate * _u[_qp];
  
  // 源项残差: -test*(G ̇ - R ̇)
  // 负号是因为将源项移到方程右侧
  return -_test[_i][_qp] * (gas_generation_rate - resolution_rate);
}