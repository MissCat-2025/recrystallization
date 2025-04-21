//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "RecrystallizationTestApp.h"
#include "RecrystallizationApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
RecrystallizationTestApp::validParams()
{
  InputParameters params = RecrystallizationApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

RecrystallizationTestApp::RecrystallizationTestApp(InputParameters parameters) : MooseApp(parameters)
{
  RecrystallizationTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

RecrystallizationTestApp::~RecrystallizationTestApp() {}

void
RecrystallizationTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  RecrystallizationApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"RecrystallizationTestApp"});
    Registry::registerActionsTo(af, {"RecrystallizationTestApp"});
  }
}

void
RecrystallizationTestApp::registerApps()
{
  registerApp(RecrystallizationApp);
  registerApp(RecrystallizationTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
RecrystallizationTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  RecrystallizationTestApp::registerAll(f, af, s);
}
extern "C" void
RecrystallizationTestApp__registerApps()
{
  RecrystallizationTestApp::registerApps();
}
