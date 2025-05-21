################################################################################
## 模型step4：简化的耦合位错密度演化、气泡形成与晶粒生长
## mpirun -n 9 /home/yp/projects/recrystallization/recrystallization-opt -i step3.3.i
################################################################################

# === 重结晶模拟 ===
# 模型step3.3：简化的耦合位错密度演化、气泡形成与晶粒生长

# 长度单位转换
length_scale = 1e-6  # 转换因子：1 = 1e-6 m = 1 μm
domain_size = 1   # 模拟域大小 (μm)

# 网格参数
nx = 50
ny = 50

# 位错密度参数 (SI单位制)
shear_modulus = 36.0e9          # 剪切模量 (Pa)=kg/(m·s²)
burgers_vector_SI = 3.42e-10    # 伯格斯矢量 (m)
rho_d_0_SI = 6.0e13             # 初始位错密度 (m^-2)
lambda_SI = 2.0e-27             # 位错演化系数 (m^3)
fission_rate_SI = 5.0e20        # 裂变率 (m^-3 s^-1)

# 气泡能量和相场参数
A = 3.0e7                       # 系数 A J(m^3)
C_p = 1.5                       # 系数 C_p
C_q = 1.8                       # 系数 C_q
S_I = 1                         # 界面能系数
c_m_e = 1e-7                    # 材料中平衡气体浓度
c_b_e = 1                       # 气泡中平衡气体浓度
initial_gas_conc = 0.012        # 初始气体浓度

# 动力学参数 (SI单位制)
L_eta_SI = 1.82e-14             # 序参数迁移率 (m^3/J·s)
kappa_eta_SI = 3.75e-8          # 序参数梯度系数 (J/m)
kappa_c_SI = 2.74e-7            # 浓度梯度系数 (J/m)
M_b_SI = 2.5e-25                # 气体浓度迁移率 (m^5/J·s)
M_g_SI = 2.5e-23              # 晶界处气体浓度迁移率 (m^5/J·s)
# 气泡形成参数 (SI单位制)
bubble_resolution_rate = 2.0e-4 # 气泡解析率 (s^-1)
atomic_volume_SI = 2.0e-29      # 原子体积 (m^3)

# 单位转换为微米 (μm)
burgers_vector = '${fparse burgers_vector_SI/length_scale}'  # (μm)
rho_d_0 = '${fparse rho_d_0_SI*length_scale*length_scale}'  # (μm^-2)
lambda = '${fparse lambda_SI/(length_scale*length_scale*length_scale)}'  # (μm^3)
fission_rate = '${fparse fission_rate_SI*length_scale*length_scale*length_scale}'  # (μm^-3 s^-1)
L_eta = '${fparse L_eta_SI/(length_scale*length_scale*length_scale)}'  # (μm^3/J·s)
kappa_eta = '${fparse kappa_eta_SI*length_scale}'  # (J/μm)
kappa_c = '${fparse kappa_c_SI*length_scale}'  # (J/μm)
M_b = '${fparse M_b_SI/(length_scale*length_scale*length_scale*length_scale*length_scale)}'  # (μm^5/J·s)
atomic_volume = '${fparse atomic_volume_SI/(length_scale*length_scale*length_scale)}'  # (μm^3)


M_g = '${fparse M_g_SI/(length_scale*length_scale*length_scale*length_scale*length_scale)}'  # (μm^5/J·s)

[Mesh]
  type        = GeneratedMesh
  dim         = 2
  nx          = ${nx}
  ny          = ${ny}
  xmin        = 0
  xmax        = ${domain_size}  # μm
  ymin        = 0
  ymax        = ${domain_size}  # μm
[]
  
[GlobalParams]
  outputs     = exodus
  derivative_order = 2
[]
  
#------------------------------------------------------------------------------#
# 1. 辅助变量：位错密度, 气泡成核概率, 局部能量
#------------------------------------------------------------------------------#
[AuxVariables]
  [./bnds]   # 晶界指示
    order     = CONSTANT
    family    = MONOMIAL
  [../]
[]
  
[AuxKernels]
  [./bnds_auxiliary]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'timestep_end'
    v = 'eta1 bubble_eta1'
  [../]
[]
  
#------------------------------------------------------------------------------#
# 2. 主变量：晶粒相场、气泡相场、气体浓度
#------------------------------------------------------------------------------#
[Variables]
  # 晶粒相场变量
  [./eta1]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = RandomIC
      seed     = 1
    [../]
  [../]
  
  # 气泡相场变量
  [./bubble_eta1]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = ConstantIC
      value    = 0.05    # 初始无气泡
    [../]
  [../]
  
  # 气体浓度变量 c_g
  [./c_g]
    order      = FIRST
    family     = LAGRANGE
    [./InitialCondition]
      type     = ConstantIC
      value    = ${initial_gas_conc}  # 初始气体浓度很低
    [../]
  [../]
[]
  
#------------------------------------------------------------------------------#
# 3. Kernels：Allen-Cahn + Cahn-Hilliard
#------------------------------------------------------------------------------#
[Kernels]
  [source_sink]
    type = ADGasSourceSink
    variable = c_g
    fission_rate = ${fission_rate}     # μm^(-3) s^(-1)
    atomic_volume = ${atomic_volume}   # μm^3
    resolution_rate = ${bubble_resolution_rate}  # s^(-1)
    bubble_etas = 'bubble_eta1'
    seed = 1234
    update_seed_per_timestep = true
  []
  #-------------------#
  # 晶粒相场变量方程 #
  #-------------------#
  # eta1
  [./deta1dt]
    type      = TimeDerivative
    variable  = eta1
  [../]
  [./ACBulk1]
    type      = AllenCahn
    variable  = eta1
    f_name    = F_total
    mob_name  = L_eta
  [../]
  
  [./ACInterface1]
    type      = ACInterface
    variable  = eta1
    mob_name  = L_eta
    kappa_name = kappa_eta
  [../]

  #-------------------#
  # 气泡相场变量方程 #
  #-------------------#
  # bubble_eta1
  [./dbubble_etadt1]
    type      = TimeDerivative
    variable  = bubble_eta1
  [../]
  
  [./ACBulk_bubble1]
    type      = AllenCahn
    variable  = bubble_eta1
    f_name    = F_total
    mob_name  = L_eta
  [../]
  
  [./ACInterface_bubble1]
    type      = ACInterface
    variable  = bubble_eta1
    mob_name  = L_eta
    kappa_name = kappa_eta
  [../]

  #-------------------#
  # 气体浓度变量方程 #
  #-------------------#
  # c_g
  [./dc_gdt]
    type      = TimeDerivative
    variable  = c_g
  [../]
  
  [./CHBulk]
    type      = CahnHilliard
    variable  = c_g
    f_name    = F_total
    mob_name  = M_c_g
  [../]
  
  [./CHInterface]
    type      = CHInterface
    variable  = c_g
    mob_name  = M_c_g
    kappa_name = kappa_c
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
  # 物理常数 - 调整参数提高稳定性
  # 气体浓度迁移率 M(η) = M_bulk * h(η) + M_gb * [1-h(η)]
  [./consts]
    type         = GenericConstantMaterial
    prop_names   = 'L_eta kappa_eta kappa_c M_b M_g' 
    prop_values  = '${L_eta} ${kappa_eta} ${kappa_c} ${M_b} ${M_g}'
  [../]

  [disloc_energy]
    type = ADDislocEnergyDensity
    # 材料参数
    shear_modulus = ${shear_modulus}  # Pa
    burgers_vector = ${burgers_vector}  # μm
    rho_d_0 = ${rho_d_0}  # μm^-2
    lambda = ${lambda}
    fission_rate = ${fission_rate}  # μm^-3 s^-1
    # 相场变量
    etas = 'eta1'
    # 输出属性名称
    dislocation_density_name = dislocation_density
    f_name = f_stored
  []
  
  [bulk_free_energy]
    type = ADBulkFreeEnergyDensity3
    # 材料相关参数
    A = ${A}
    C_p = ${C_p}
    C_q = ${C_q}
    S_I = ${S_I}
    c_m_e = ${c_m_e}
    c_b_e = ${c_b_e}
    # 变量
    c_g = c_g  # 气体浓度变量
    grain_etas = 'eta1'  # 晶粒序参数变量
    bubble_etas = 'bubble_eta1'  # 气泡序参数变量
    
    # 输出属性名称
    f_name = f_bulk
  []
  
  # 晶粒切换函数 h_i
  [./switching1]
    type         = SwitchingFunctionMaterial
    eta          = eta1
    h_order      = HIGH
  [../]
  
  # 总自由能 F_total = F_i + F_interaction + F_gas + 0.25
  [./total_free_energy]
    type             = DerivativeParsedMaterial
    property_name    = F_total
    coupled_variables= 'eta1 bubble_eta1 c_g'
    expression       = 'f_bulk + f_stored'
    material_property_names = 'f_bulk f_stored'
    outputs          = exodus
  [../]

  [gas_mobility]
    type = DerivativeParsedMaterial
    property_name = M_c_g
    material_property_names = 'h(eta1) M_b M_g'
    expression = 'M_b*h + M_g*(1-h)'
    outputs = exodus
  []
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
    dt = 0.01  # 使用非常小的初始时间步长
  [../]
  nl_rel_tol      = 1e-7
  nl_abs_tol      = 1e-9
  nl_max_its      = 25
  l_max_its       = 100
  l_tol           = 1e-6
  line_search     = 'bt'  # 使用回溯线搜索
  [./Adaptivity]
    initial_adaptivity = 2
    refine_fraction = 0.3
    coarsen_fraction = 0.05
  [../]
  automatic_scaling = true
  compute_scaling_once = false
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 200'
[]
  
[Outputs]
  [./exodus]
    type = Exodus
    time_step_interval = 1
    # show = 'eta1 bubble_eta1 c_g bnds'
  [../]
  print_linear_residuals = true
  checkpoint = true
[]