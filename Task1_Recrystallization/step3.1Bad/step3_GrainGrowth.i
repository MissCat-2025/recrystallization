################################################################################
## 模型step2：耦合位错密度演化、气泡形成与再结晶晶粒生长
## mpirun -n 12 /home/yp/projects/recrystallization/recrystallization-opt -i step3_GrainGrowth.i
################################################################################

[Mesh]
  type        = GeneratedMesh
  dim         = 2
  nx          = 40
  ny          = 40
  xmin        = 0
  xmax        = 10
  ymin        = 0
  ymax        = 10
[]

[GlobalParams]
  outputs     = exodus
  penalty     = 500
  derivative_order = 2
[]

#------------------------------------------------------------------------------#
# 1. 辅助变量：位错密度, 气泡成核概率, 局部能量
#------------------------------------------------------------------------------#
[AuxVariables]
  [./dislocation_density]  # 位错密度 ρ_d(t)
    order      = CONSTANT
    family     = MONOMIAL
    [./InitialCondition]
      type     = ConstantIC
      value    = 1.0e12    # 初始位错密度 ρ_d(0)
    [../]
  [../]
  [./nucleation_probability]  # 气泡成核概率
    order      = CONSTANT
    family     = MONOMIAL
  [../]
  [./stored_energy]  # 位错储存弹性能
    order      = CONSTANT
    family     = MONOMIAL
  [../]
  [./local_energy]  # 局部自由能
    order      = CONSTANT
    family     = MONOMIAL
  [../]
  [./cross_energy]  # 梯度交叉项
    order      = CONSTANT
    family     = MONOMIAL
  [../]
  [./bnds]   # 晶界指示
    order     = CONSTANT
    family    = MONOMIAL
  [../]
  [./bubble_sum_squared]   # 气泡相场平方和
    order     = CONSTANT
    family    = MONOMIAL
  [../]
  []

[Functions]
  [./dislocation_growth]
    type = ParsedFunction
    # 直接计算出表达式值，避免使用未定义的变量u
    expression = 'exp(4.0e3 * 5.0e20 * 0.05)'  # 常数表达式，计算增长因子
  [../]
[]

[AuxKernels]
  # 修改位错密度演化计算方式
  [./dislocation_evolution]
    type      = ParsedAux
    variable  = dislocation_density
    coupled_variables = 'dislocation_density'
    expression = 'dislocation_density * exp(4.0e3 * 5.0e20 * 0.05)'  # 直接计算，无需使用函数
    execute_on = 'timestep_begin'
  [../]
  
  # 计算储存弹性能 f_stored(η_i) = 1/2·η_i^2·G·b_v^2·ρ_d
  [./stored_elastic_energy]
    type      = ParsedAux
    variable  = stored_energy
    coupled_variables = 'dislocation_density eta1 eta2 eta3 eta4 eta5 eta6'
    expression = '0.5 * 36.0e9 * (3.42e-10)^2 * dislocation_density * (eta1^2 + eta2^2 + eta3^2 + eta4^2 + eta5^2 + eta6^2)'
  [../]
  
  # 计算局部自由能
  [./local_free_energy]
    type               = ParsedAux  # 改为ParsedAux计算自由能
    variable           = local_energy
    coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    expression         = 'eta1^2 + eta2^2 + eta3^2 + eta4^2 + eta5^2 + eta6^2 + bubble_eta1^2 + bubble_eta2^2 + bubble_eta3^2 + bubble_eta4^2'  # 简化的自由能表达式
  [../]
  
  # 计算梯度交叉项
  [./cross_terms]
    type               = CrossTermGradientFreeEnergy
    variable           = cross_energy
    interfacial_vars   = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    # 10x10矩阵的所有组合
    kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble
                          kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble
                          kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble
                          kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble
                          kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble
                          kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble
                          kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble
                          kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble
                          kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble
                          kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
  [../]
  
  # 计算晶界
  [./bnds_auxiliary]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'timestep_end'
    v = 'eta1 eta2 eta3 eta4 eta5 eta6'
  [../]
  
  # 计算气泡成核概率 P(t) = 1 - exp(-J·Δt)
  [./bubble_nucleation_probability]
    type = ParsedAux
    variable = nucleation_probability
    coupled_variables = 'bnds c_g'
    expression = 'if(bnds > 0.01, 1 - exp(-5.0e-3 * exp(-1.0e-6/(c_g - 1.0e-7))), 0)'
  [../]
  [./bubble_sum_squared_aux]
    type = ParsedAux
    variable = bubble_sum_squared
    coupled_variables = 'bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    expression = 'bubble_eta1^2 + bubble_eta2^2 + bubble_eta3^2 + bubble_eta4^2'
    execute_on = 'timestep_begin timestep_end'
  [../]
  []

#------------------------------------------------------------------------------#
# 2. 主变量：晶粒相场、气泡相场、气体浓度
#------------------------------------------------------------------------------#
[Variables]
  # 晶粒相场变量 (i = 1, 2, ..., p)
  [./eta1]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = RandomIC
      seed     = 1
    [../]
  [../]
  [./eta2]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = RandomIC
      seed     = 2
    [../]
  [../]
  [./eta3]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = RandomIC
      seed     = 3
    [../]
  [../]
  [./eta4]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = RandomIC
      seed     = 4
    [../]
  [../]
  [./eta5]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = RandomIC
      seed     = 5
    [../]
  [../]
  [./eta6]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = RandomIC
      seed     = 6
    [../]
  [../]
  
  # 气泡相场变量 (j = p+1, ..., q)
  [./bubble_eta1]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = ConstantIC
      value    = 0.0    # 初始无气泡
    [../]
  [../]
  [./bubble_eta2]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = ConstantIC
      value    = 0.0    # 初始无气泡
    [../]
  [../]
  [./bubble_eta3]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = ConstantIC
      value    = 0.0    # 初始无气泡
    [../]
  [../]
  [./bubble_eta4]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = ConstantIC
      value    = 0.0    # 初始无气泡
    [../]
  [../]
  
  # 气体浓度变量 c_g
  [./c_g]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = ConstantIC
      value    = 0.001  # 初始气体浓度很低
    [../]
  [../]
  []

#------------------------------------------------------------------------------#
# 3. Kernels：Allen-Cahn + ACMultiInterface + Cahn-Hilliard
#------------------------------------------------------------------------------#
[Kernels]
  #-------------------#
  # 晶粒相场变量方程 #
  #-------------------#
  # eta1
  [./deta1dt]
    type      = TimeDerivative
    variable  = eta1
  [../]
  [./ACBulk1]
    type               = AllenCahn
    variable           = eta1
    f_name             = F_total
    mob_name           = L_eta
  [../]
  [./ACInterface1]
    type               = ACMultiInterface
    variable           = eta1
    etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    mob_name           = L_eta
    kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
  [../]
  [./penalty1]
    type      = SwitchingFunctionPenalty
    variable  = eta1
    etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    h_names   = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4'
  [../]
  [./elastic_energy1]
    type      = AllenCahn
    variable  = eta1
    f_name    = F_stored
    mob_name  = L_eta
    coupled_variables = 'eta2 eta3 eta4 eta5 eta6'
  [../]

  # eta2
  [./deta2dt]
    type      = TimeDerivative
    variable  = eta2
  [../]
  [./ACBulk2]
    type               = AllenCahn
    variable           = eta2
    f_name             = F_total
    mob_name           = L_eta
  [../]
  [./ACInterface2]
    type               = ACMultiInterface
    variable           = eta2
    etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    mob_name           = L_eta
    kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
  [../]
  [./penalty2]
    type      = SwitchingFunctionPenalty
    variable  = eta2
    etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    h_names   = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4'
  [../]
  [./elastic_energy2]
    type      = AllenCahn
    variable  = eta2
    f_name    = F_stored
    mob_name  = L_eta
    coupled_variables = 'eta1 eta3 eta4 eta5 eta6'
  [../]

  # eta3
  [./deta3dt]
    type      = TimeDerivative
    variable  = eta3
  [../]
  [./ACBulk3]
    type               = AllenCahn
    variable           = eta3
    f_name             = F_total
    mob_name           = L_eta
  [../]
  [./ACInterface3]
    type               = ACMultiInterface
    variable           = eta3
    etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    mob_name           = L_eta
    kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
  [../]
  [./penalty3]
    type      = SwitchingFunctionPenalty
    variable  = eta3
    etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    h_names   = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4'
  [../]
  [./elastic_energy3]
    type      = AllenCahn
    variable  = eta3
    f_name    = F_stored
    mob_name  = L_eta
    coupled_variables = 'eta1 eta2 eta4 eta5 eta6'
  [../]

  # eta4
  [./deta4dt]
    type      = TimeDerivative
    variable  = eta4
  [../]
  [./ACBulk4]
    type               = AllenCahn
    variable           = eta4
    f_name             = F_total
    mob_name           = L_eta
  [../]
  [./ACInterface4]
    type               = ACMultiInterface
    variable           = eta4
    etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    mob_name           = L_eta
    kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
  [../]
  [./penalty4]
    type      = SwitchingFunctionPenalty
    variable  = eta4
    etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    h_names   = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4'
  [../]
  [./elastic_energy4]
    type      = AllenCahn
    variable  = eta4
    f_name    = F_stored
    mob_name  = L_eta
    coupled_variables = 'eta1 eta2 eta3 eta5 eta6'
  [../]

  # eta5
  [./deta5dt]
    type      = TimeDerivative
    variable  = eta5
  [../]
  [./ACBulk5]
    type               = AllenCahn
    variable           = eta5
    f_name             = F_total
    mob_name           = L_eta
  [../]
  [./ACInterface5]
    type               = ACMultiInterface
    variable           = eta5
    etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    mob_name           = L_eta
    kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
  [../]
  [./penalty5]
    type      = SwitchingFunctionPenalty
    variable  = eta5
    etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    h_names   = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4'
  [../]
  [./elastic_energy5]
    type      = AllenCahn
    variable  = eta5
    f_name    = F_stored
    mob_name  = L_eta
    coupled_variables = 'eta1 eta2 eta3 eta4 eta6'
  [../]

  # eta6
  [./deta6dt]
    type      = TimeDerivative
    variable  = eta6
  [../]
  [./ACBulk6]
    type               = AllenCahn
    variable           = eta6
    f_name             = F_total
    mob_name           = L_eta
  [../]
  [./ACInterface6]
    type               = ACMultiInterface
    variable           = eta6
    etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    mob_name           = L_eta
    kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
  [../]
  [./penalty6]
    type      = SwitchingFunctionPenalty
    variable  = eta6
    etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    h_names   = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4'
  [../]
  [./elastic_energy6]
    type      = AllenCahn
    variable  = eta6
    f_name    = F_stored
    mob_name  = L_eta
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5'
  [../]

  #-------------------#
  # 气泡相场变量方程 #
  #-------------------#
  [./dbubble_etadt1]
    type      = TimeDerivative
    variable  = bubble_eta1
  [../]
  [./ACBulk_bubble1]
    type               = AllenCahn
    variable           = bubble_eta1
    f_name             = F_total
    mob_name           = L_bubble
  [../]
  [./ACInterface_bubble1]
    type               = ACMultiInterface
    variable           = bubble_eta1
    etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    mob_name           = L_bubble
    kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
  [../]
  [./penalty_bubble1]
    type      = SwitchingFunctionPenalty
    variable  = bubble_eta1
    etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    h_names   = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4'
  [../]
  [./bubble_nucleation1]
    type      = LangevinNoise
    variable  = bubble_eta1
    amplitude = 5.0  # 直接使用数值而非nucleation_factor
  [../]
  
  #气泡2
    [./dbubble_etadt2]
      type      = TimeDerivative
      variable  = bubble_eta2
    [../]
    [./ACBulk_bubble2]
      type               = AllenCahn
      variable           = bubble_eta2
      f_name             = F_total
      mob_name           = L_bubble
    [../]
    [./ACInterface_bubble2]
      type               = ACMultiInterface
      variable           = bubble_eta2
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
      mob_name           = L_bubble
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
    [../]
    [./penalty_bubble2]
      type      = SwitchingFunctionPenalty
      variable  = bubble_eta2
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
      h_names   = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4'
    [../]
    [./bubble_nucleation2]
      type      = LangevinNoise
      variable  = bubble_eta2
      amplitude = 5.0  # 直接使用数值而非nucleation_factor
    [../]

      #气泡3
    [./dbubble_etadt3]
      type      = TimeDerivative
      variable  = bubble_eta3
    [../]
    [./ACBulk_bubble3]
      type               = AllenCahn
      variable           = bubble_eta3
      f_name             = F_total
      mob_name           = L_bubble
    [../]
    [./ACInterface_bubble3]
      type               = ACMultiInterface
      variable           = bubble_eta3
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
      mob_name           = L_bubble
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
    [../]
    [./penalty_bubble3]
      type      = SwitchingFunctionPenalty
      variable  = bubble_eta3
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
      h_names   = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4'
    [../]
    [./bubble_nucleation3]
      type      = LangevinNoise
      variable  = bubble_eta3
      amplitude = 5.0  # 直接使用数值而非nucleation_factor
    [../]

      #气泡4
    [./dbubble_etadt4]
      type      = TimeDerivative
      variable  = bubble_eta4
    [../]
    [./ACBulk_bubble4]
      type               = AllenCahn
      variable           = bubble_eta4
      f_name             = F_total
      mob_name           = L_bubble
    [../]
    [./ACInterface_bubble4]
      type               = ACMultiInterface
      variable           = bubble_eta4
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
      mob_name           = L_bubble
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_bubble kappa_bubble kappa_bubble kappa_bubble'
    [../]
    [./penalty_bubble4]
      type      = SwitchingFunctionPenalty
      variable  = bubble_eta4
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
      h_names   = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4'
    [../]
    [./bubble_nucleation4]
      type      = LangevinNoise
      variable  = bubble_eta4
      amplitude = 5.0  # 直接使用数值而非nucleation_factor
    [../]
  #-------------------#
  # 气体浓度变量方程 #
  #-------------------#
  [./dc_gdt]
    type      = TimeDerivative
    variable  = c_g
  [../]
  [./CHBulk]
    type      = CahnHilliard
    variable  = c_g
    f_name    = F_gas
    mob_name  = M_c_g
  [../]
  [./CHInterface]
    type      = CHInterface
    variable  = c_g
    mob_name  = M_c_g
    kappa_name = kappa_c
  [../]
  [./gas_generation]
    type      = BodyForce
    variable  = c_g
    value     = 0.012  # 气体产生率 G_dot
  [../]
  [./gas_resolution]
    type      = CoefReaction
    variable  = c_g
    # reaction_rate = -2.0e-4  # 直接使用数值代替材料属性
    # v = bubble_sum_squared
  [../]
  [./bubble1_concentration_coupling]
    type      = CoupledSwitchingTimeDerivative
    variable  = bubble_eta1
    v         = c_g
    # eta       = bubble_eta1
    # h_name    = h_bubble1
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6'
    Fj_names  = 'F1 F2 F3 F4 F5 F6'
    hj_names  = 'h1 h2 h3 h4 h5 h6'
  [../]
  [./bubble1_grain_interaction]
    type      = AllenCahn
    variable  = bubble_eta1
    f_name    = F_grain_bubble_interaction
    mob_name  = L_bubble
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta2 bubble_eta3 bubble_eta4'
  [../]
  
  [./bubble2_concentration_coupling]
    type      = CoupledSwitchingTimeDerivative
    variable  = bubble_eta2
    v         = c_g
    # eta       = bubble_eta2
    # h_name    = h_bubble2
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6'
    Fj_names  = 'F1 F2 F3 F4 F5 F6'
    hj_names  = 'h1 h2 h3 h4 h5 h6'
  [../]
  [./bubble2_grain_interaction]
    type      = AllenCahn
    variable  = bubble_eta2
    f_name    = F_grain_bubble_interaction
    mob_name  = L_bubble
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta3 bubble_eta4'
  [../]
  
  [./bubble3_concentration_coupling]
    type      = CoupledSwitchingTimeDerivative
    variable  = bubble_eta3
    v         = c_g
      # eta       = bubble_eta3
    # h_name    = h_bubble3
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6'
    Fj_names  = 'F1 F2 F3 F4 F5 F6'
    hj_names  = 'h1 h2 h3 h4 h5 h6'
  [../]
  
  [./bubble3_grain_interaction]
    type      = AllenCahn
    variable  = bubble_eta3
    f_name    = F_grain_bubble_interaction
    mob_name  = L_bubble
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta4'
  [../]
  
  [./bubble4_concentration_coupling]
    type      = CoupledSwitchingTimeDerivative
    variable  = bubble_eta4
    v         = c_g
    # eta       = bubble_eta4
    # h_name    = h_bubble4
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6'
    Fj_names  = 'F1 F2 F3 F4 F5 F6'
    hj_names  = 'h1 h2 h3 h4 h5 h6'
  [../]
  
  [./bubble4_grain_interaction]
    type      = AllenCahn
    variable  = bubble_eta4
    f_name    = F_grain_bubble_interaction
    mob_name  = L_bubble
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3'
  [../]
  []

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y'
    [../]
  [../]
  []

#------------------------------------------------------------------------------#
# 4. Materials：常数、自由能、动力学参数
#------------------------------------------------------------------------------#
[Materials]
  # 物理常数
  [./consts]
    type         = GenericConstantMaterial
    prop_names   = 'A Cp Cq L_eta L_bubble M_bulk M_gb
                    kappa_eta kappa_bubble kappa_c
                    G b_v lambda fission_rate
                    k1 k2 c_m_e c_b_e S_I
                    gas_gen_rate resolution_rate
                    nucleation_factor dt T'
    prop_values  = '3.0e7 1.5 1.8 1.82e-14 1.82e-14 2.5e-25 2.5e-23
                    3.75e-8 3.75e-8 2.74e-7
                    36.0e9 3.42e-10 4.0e3 5.0e20
                    5.0e-3 1.0e-6 1.0e-7 1.0 1.0
                    0.012 2.0e-4
                    5.0 0.05 473'  # 更新参数值
  [../]
  
  # 添加一个注释材料，便于跟踪参数值
  [./params_info]
    type         = ParsedMaterial
    property_name = params_info
    material_property_names = 'A Cp Cq L_eta M_bulk M_gb kappa_eta kappa_c G b_v fission_rate'
    expression   = '0' # 只是为了记录参数，不参与计算
    outputs      = exodus
  [../]
  
  # 气体浓度迁移率 M(η) = M_bulk * h(η) + M_gb * [1-h(η)]
  [./gas_mobility]
    type         = DerivativeParsedMaterial
    property_name = M_c_g
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6'
    expression   = 'M_bulk * (sum_eta^3 * (6*sum_eta^2 - 15*sum_eta + 10)) + 
                    M_gb * (1 - (sum_eta^3 * (6*sum_eta^2 - 15*sum_eta + 10)))'
    material_property_names = 'M_bulk M_gb sum_eta'
    constant_on  = 'ELEMENT'
    outputs      = exodus
  [../]
  
  # 辅助材料属性，计算总相场值
  [./sum_eta_aux]
    type = ParsedMaterial
    property_name = sum_eta
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6'
    expression = 'eta1 + eta2 + eta3 + eta4 + eta5 + eta6'
    outputs = exodus
  [../]
  
  # 阻挡函数 ∑ h_i = 1
  [./barrier]
    type         = MultiBarrierFunctionMaterial
    etas         = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
  [../]
  
  # 晶粒切换函数 h_i
  [./switching1]
    type         = SwitchingFunctionMaterial
    function_name= h1
    eta          = eta1
    h_order      = SIMPLE
  [../]
  [./switching2]
    type         = SwitchingFunctionMaterial
    function_name= h2
    eta          = eta2
    h_order      = SIMPLE
  [../]
  [./switching3]
    type         = SwitchingFunctionMaterial
    function_name= h3
    eta          = eta3
    h_order      = SIMPLE
  [../]
  [./switching4]
    type         = SwitchingFunctionMaterial
    function_name= h4
    eta          = eta4
    h_order      = SIMPLE
  [../]
  [./switching5]
    type         = SwitchingFunctionMaterial
    function_name= h5
    eta          = eta5
    h_order      = SIMPLE
  [../]
  [./switching6]
    type         = SwitchingFunctionMaterial
    function_name= h6
    eta          = eta6
    h_order      = SIMPLE
  [../]
  
  # 气泡切换函数
  [./switching_bubble1]
    type         = SwitchingFunctionMaterial
    function_name= h_bubble1
    eta          = bubble_eta1
    h_order      = SIMPLE
  [../]
  
  [./switching_bubble2]
    type         = SwitchingFunctionMaterial
    function_name= h_bubble2
    eta          = bubble_eta2
    h_order      = SIMPLE
  [../]
  
  [./switching_bubble3]
    type         = SwitchingFunctionMaterial
    function_name= h_bubble3
    eta          = bubble_eta3
    h_order      = SIMPLE
  [../]
  
  [./switching_bubble4]
    type         = SwitchingFunctionMaterial
    function_name= h_bubble4
    eta          = bubble_eta4
    h_order      = SIMPLE
  [../]
  
  # 体积分数计算
  [./theta_m]
    type         = DerivativeParsedMaterial
    property_name = theta_m
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    expression   = '(eta1^2 + eta2^2 + eta3^2 + eta4^2 + eta5^2 + eta6^2) / 
                    (eta1^2 + eta2^2 + eta3^2 + eta4^2 + eta5^2 + eta6^2 + bubble_eta1^2 + bubble_eta2^2 + bubble_eta3^2 + bubble_eta4^2)'
    outputs      = exodus
  [../]
  
  [./theta_b]
    type         = DerivativeParsedMaterial
    property_name = theta_b
    coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    expression   = '1 - theta_m'
    material_property_names = 'theta_m'
    outputs      = exodus
  [../]
  
  # 气体自由能 f_g^m = S_I * (c_g - c_m_e)^2, f_g^b = S_I * (c_g - c_b_e)^2
  [./gas_free_energy_m]
    type         = DerivativeParsedMaterial
    property_name = f_g_m
    coupled_variables = 'c_g'
    expression   = 'S_I * (c_g - c_m_e)^2'
    material_property_names = 'S_I c_m_e'
  [../]
  
  [./gas_free_energy_b]
    type         = DerivativeParsedMaterial
    property_name = f_g_b
    coupled_variables = 'c_g'
    expression   = 'S_I * (c_g - c_b_e)^2'
    material_property_names = 'S_I c_b_e'
  [../]
  
  # 位错储存能 f_stored(η_i) = 1/2·η_i^2·G·b_v^2·ρ_d
  [./stored_elastic_energy]
    type         = DerivativeParsedMaterial
    property_name = F_stored
    coupled_variables = 'dislocation_density eta1 eta2 eta3 eta4 eta5 eta6'
    expression   = '0.5 * G * b_v^2 * dislocation_density * (eta1^2 + eta2^2 + eta3^2 + eta4^2 + eta5^2 + eta6^2)'
    material_property_names = 'G b_v'
    outputs      = exodus
  [../]
  
  # 各相自由能 F_i：简单双井势
  [./phase_free_energy_1]
    type             = DerivativeParsedMaterial
    property_name    = F1
    coupled_variables= 'eta1'
    expression       = 'A*(-0.5*eta1^2 + 0.25*eta1^4)'
    material_property_names = 'A'
  [../]
  [./phase_free_energy_2]
    type             = DerivativeParsedMaterial
    property_name    = F2
    coupled_variables= 'eta2'
    expression       = 'A*(-0.5*eta2^2 + 0.25*eta2^4)'
    material_property_names = 'A'
  [../]
  [./phase_free_energy_3]
    type             = DerivativeParsedMaterial
    property_name    = F3
    coupled_variables= 'eta3'
    expression       = 'A*(-0.5*eta3^2 + 0.25*eta3^4)'
    material_property_names = 'A'
  [../]
  [./phase_free_energy_4]
    type             = DerivativeParsedMaterial
    property_name    = F4
    coupled_variables= 'eta4'
    expression       = 'A*(-0.5*eta4^2 + 0.25*eta4^4)'
    material_property_names = 'A'
  [../]
  [./phase_free_energy_5]
    type             = DerivativeParsedMaterial
    property_name    = F5
    coupled_variables= 'eta5'
    expression       = 'A*(-0.5*eta5^2 + 0.25*eta5^4)'
    material_property_names = 'A'
  [../]
  [./phase_free_energy_6]
    type             = DerivativeParsedMaterial
    property_name    = F6
    coupled_variables= 'eta6'
    expression       = 'A*(-0.5*eta6^2 + 0.25*eta6^4)'
    material_property_names = 'A'
  [../]
  
  # 气泡自由能
  [./bubble_free_energy1]
    type             = DerivativeParsedMaterial
    property_name    = F_bubble1
    coupled_variables= 'bubble_eta1'
    expression       = 'A*(-0.5*bubble_eta1^2 + 0.25*bubble_eta1^4)'
    material_property_names = 'A'
  [../]
  
  [./bubble_free_energy2]
    type             = DerivativeParsedMaterial
    property_name    = F_bubble2
    coupled_variables= 'bubble_eta2'
    expression       = 'A*(-0.5*bubble_eta2^2 + 0.25*bubble_eta2^4)'
    material_property_names = 'A'
  [../]
  
  [./bubble_free_energy3]
    type             = DerivativeParsedMaterial
    property_name    = F_bubble3
    coupled_variables= 'bubble_eta3'
    expression       = 'A*(-0.5*bubble_eta3^2 + 0.25*bubble_eta3^4)'
    material_property_names = 'A'
  [../]
  
  [./bubble_free_energy4]
    type             = DerivativeParsedMaterial
    property_name    = F_bubble4
    coupled_variables= 'bubble_eta4'
    expression       = 'A*(-0.5*bubble_eta4^2 + 0.25*bubble_eta4^4)'
    material_property_names = 'A'
  [../]
  
  # 晶粒-晶粒相互作用 Cp*sum(i,j!=i) η_i^2 * η_j^2
  [./grain_interaction]
    type             = DerivativeParsedMaterial
    property_name    = F_grain_interaction
    coupled_variables= 'eta1 eta2 eta3 eta4 eta5 eta6'
    expression       = 'Cp * (
                        eta1^2 * (eta2^2 + eta3^2 + eta4^2 + eta5^2 + eta6^2) +
                        eta2^2 * (eta3^2 + eta4^2 + eta5^2 + eta6^2) +
                        eta3^2 * (eta4^2 + eta5^2 + eta6^2) +
                        eta4^2 * (eta5^2 + eta6^2) +
                        eta5^2 * eta6^2)'
    material_property_names = 'Cp'
  [../]
  
  # 晶粒-气泡相互作用 Cq*sum(i,j) η_i^2 * η_j^2
  [./grain_bubble_interaction]
    type             = DerivativeParsedMaterial
    property_name    = F_grain_bubble_interaction
    coupled_variables= 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4'
    expression       = 'Cq * (bubble_eta1^2 + bubble_eta2^2 + bubble_eta3^2 + bubble_eta4^2) * (eta1^2 + eta2^2 + eta3^2 + eta4^2 + eta5^2 + eta6^2)'
    material_property_names = 'Cq'
  [../]
  
  # 气体化学自由能 f_g^m * θ_m + f_g^b * θ_b
  [./gas_chemical_energy]
    type             = DerivativeParsedMaterial
    property_name    = F_gas
    coupled_variables= 'c_g'
    expression       = 'f_g_m * theta_m + f_g_b * theta_b'
    material_property_names = 'f_g_m f_g_b theta_m theta_b'
    outputs          = exodus
  [../]
  
  # 总自由能 F_total = F_i + F_interaction + F_gas + 0.25
  [./total_free_energy]
    type             = DerivativeParsedMaterial
    property_name    = F_total
    coupled_variables= 'eta1 eta2 eta3 eta4 eta5 eta6 bubble_eta1 bubble_eta2 bubble_eta3 bubble_eta4 c_g'
    expression       = 'h1*F1 + h2*F2 + h3*F3 + h4*F4 + h5*F5 + h6*F6 + 
                        h_bubble1*F_bubble1 + h_bubble2*F_bubble2 + h_bubble3*F_bubble3 + h_bubble4*F_bubble4 + 
                        F_grain_interaction + F_grain_bubble_interaction + 
                        F_gas + 0.25'
    material_property_names = 'h1 h2 h3 h4 h5 h6 h_bubble1 h_bubble2 h_bubble3 h_bubble4
                              F1 F2 F3 F4 F5 F6 F_bubble1 F_bubble2 F_bubble3 F_bubble4
                              F_grain_interaction F_grain_bubble_interaction
                              F_gas'
    outputs          = exodus
  [../]
  []

[Postprocessors]
  [./total_energy]
    type      = ElementIntegralVariablePostprocessor
    variable  = local_energy
  [../]
  [./bubble_volume_fraction]
    type      = ElementAverageValue
    variable  = bubble_eta1
  [../]
  [./average_gas_concentration]
    type      = ElementAverageValue
    variable  = c_g
  [../]
  [./average_dislocation_density]
    type      = ElementAverageValue
    variable  = dislocation_density
  [../]
  []

[Preconditioning]
  [./SMP]
    type      = SMP
    full      = true
    petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
    petsc_options_value = 'asm      lu           2               200'
  [../]
  []

[Executioner]
  type            = Transient
  scheme          = bdf2
  solve_type      = PJFNK
  start_time      = 0.0
  end_time        = 100.0
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 0.0001  # 减小初始时间步长以提高稳定性
    # cutback_factor = 0.5
    # growth_factor = 1.2
    # optimal_iterations = 8
    # timestep_limiting_function = 1.0e-12
  [../]
  nl_rel_tol      = 1e-6
  nl_abs_tol      = 1e-8
  nl_max_its      = 20  # 增加最大迭代次数
  l_max_its       = 50  # 增加线性求解器最大迭代次数
  l_tol           = 1e-5
  line_search     = 'bt'  # 使用回溯线搜索
  automatic_scaling = true
  []

[Outputs]
  [./exodus]
    type = Exodus
    time_step_interval = 10
    show = 'eta1 bubble_eta1 c_g dislocation_density stored_energy bnds bubble_sum_squared'
  [../]
  print_linear_residuals = false
  checkpoint = true
  []