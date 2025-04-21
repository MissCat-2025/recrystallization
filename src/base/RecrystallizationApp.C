#include "RecrystallizationApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
RecrystallizationApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

RecrystallizationApp::RecrystallizationApp(InputParameters parameters) : MooseApp(parameters)
{
  RecrystallizationApp::registerAll(_factory, _action_factory, _syntax);
}

RecrystallizationApp::~RecrystallizationApp() {}

void
RecrystallizationApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<RecrystallizationApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"RecrystallizationApp"});
  Registry::registerActionsTo(af, {"RecrystallizationApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
RecrystallizationApp::registerApps()
{
  registerApp(RecrystallizationApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
RecrystallizationApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  RecrystallizationApp::registerAll(f, af, s);
}
extern "C" void
RecrystallizationApp__registerApps()
{
  RecrystallizationApp::registerApps();
}
