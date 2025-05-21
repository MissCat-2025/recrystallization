# This simulation predicts GB migration of a 2D copper polycrystal with 100 grains represented with 8 order parameters
# Mesh adaptivity and time step adaptivity are used
# An AuxVariable is used to calculate the grain boundary locations
# Postprocessors are used to record time step and the number of grains
# mpirun -n 12 /home/yp/projects/recrystallization/recrystallization-opt -i grain_growth_2D_graintracker.i 
[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 44 # Number of elements in the x-direction
  ny = 44 # Number of elements in the y-direction
  xmax = 1000 # maximum x-coordinate of the mesh
  ymax = 1000 # maximum y-coordinate of the mesh
  elem_type = QUAD4 # Type of elements used in the mesh
  uniform_refine = 2 # Initial uniform refinement of the mesh
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 8 # Number of order parameters used
  var_name_base = gr # Base name of grains
[]

[Modules]
  [PhaseField]
    [GrainGrowth]
    []
  []
[]

[UserObjects]
  [voronoi]
    type = PolycrystalVoronoi
    grain_num = 20 # Number of grains
    rand_seed = 10
    int_width = 7
  []
  [grain_tracker]
    type = GrainTracker
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
[]

[Variables]
  [./c]
    order = FIRST
    family = LAGRANGE
    #
    # We set up a smooth cradial concentrtaion gradient
    # The concentration will quickly change to adapt to the preset order
    # parameters eta1, eta2, and eta3
    #
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 0.0
      y1 = 0.0
      radius = 100.0
      invalue = 1.0
      outvalue = 0.01
      int_width = 10.0
    [../]
  [../]
  []

[AuxVariables]
  # Dependent variables
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  [bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds #计算晶界(Grain Boundary)位置，通过计算相邻晶粒之间的界面能量梯度。
    execute_on = 'initial timestep_end'
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains #为每个独立的晶粒区域分配一个唯一的ID号。使用grain_tracker用户对象识别并追踪单个晶粒。
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices #显示每个区域使用的顺序参数变量（order parameter）的索引颜色。在多相场模型中，这可以显示哪个顺序参数变量主导了特定区域。
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
  []
[]

[BCs]
  # Boundary Condition block
  [Periodic]
    [All]
      auto_direction = 'x y' # Makes problem periodic in the x and y directions
    []
  []
[]

[Materials]
  [CuGrGr]
    # Material properties
    type = GBEvolution
    T = 450 # Constant temperature of the simulation (for mobility calculation)
    wGB = 14 # Width of the diffuse GB
    GBmob0 = 2.5e-6 #m^4(Js) for copper from schonfelder1997molecular bibtex entry
    Q = 0.23 #eV for copper from schonfelder1997molecular bibtex entry
    GBenergy = 0.708 #J/m^2 from schonfelder1997molecular bibtex entry
  []
[]

[Postprocessors]
  # Scalar postprocessors
  [dt]
    # Outputs the current time step
    type = TimestepSize
  []
[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # 指定时间积分格式为二阶后向差分格式(backward differentiation formula)。

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  l_max_its = 50 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 10 # Max number of nonlinear iterations

  end_time = 4000

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 20 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # 使每个时间步的非线性求解需要大约6次迭代才能收敛。
  []

  [Adaptivity]
    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.8 # Fraction of high error that will be refined
    coarsen_fraction = 0.05 # Fraction of low error that will coarsened
    max_h_level = 2 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  []
[]

[Outputs]
  exodus = true
  file_base = 'results/grain_growth'  # 自定义文件名基础
[]

# mpirun -n 12 /home/yp/projects/recrystallization/recrystallization-opt -i grain_growth_2D_graintracker.i 