#include "matrixFreePDE.h"

using namespace dealii;

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  // Constructor
  customPDE(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
    , userInputs(_userInputs) {};

  // Function to set the initial conditions (in ICs_and_BCs.h)
  void
  setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                      [[maybe_unused]] const unsigned int index,
                      [[maybe_unused]] double            &scalar_IC,
                      [[maybe_unused]] Vector<double>    &vector_IC) override;

  // Function to set the non-uniform Dirichlet boundary conditions (in
  // ICs_and_BCs.h)
  void
  setNonUniformDirichletBCs([[maybe_unused]] const Point<dim>  &p,
                            [[maybe_unused]] const unsigned int index,
                            [[maybe_unused]] const unsigned int direction,
                            [[maybe_unused]] const double       time,
                            [[maybe_unused]] double            &scalar_BC,
                            [[maybe_unused]] Vector<double>    &vector_BC) override;

private:
#include "typeDefs.h"

  const userInputParameters<dim> userInputs;

  // Function to set the RHS of the governing equations for explicit time
  // dependent equations (in equations.h)
  void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set the RHS of the governing equations for all other equations
  // (in equations.h)
  void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set the LHS of the governing equations (in equations.h)
  void
  equationLHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

// Function to set postprocessing expressions (in postprocess.h)
#ifdef POSTPROCESS_FILE_EXISTS
  void
  postProcessedFields(
    [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
      &variable_list,
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &pp_variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;
#endif

// Function to set the nucleation probability (in nucleation.h)
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability([[maybe_unused]] variableValueContainer variable_value,
                           [[maybe_unused]] double                 dV) const override;
#endif

  // ================================================================
  // Methods specific to this subclass
  // ================================================================

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

    // physical constants
    double FarC = 96485.33289;
    double R = 8.314;

    double icorr = userInputs.get_model_constant_double("icorr");
    double W = userInputs.get_model_constant_double("W");
    double epssq = userInputs.get_model_constant_double("epssq");
    double cB_eq = userInputs.get_model_constant_double("cB_eq");
    double cB_init = userInputs.get_model_constant_double("cB_init");
    double cV_init = userInputs.get_model_constant_double("cV_init");
    double T = userInputs.get_model_constant_double("T");
    double lthresh = userInputs.get_model_constant_double("lthresh");
    double DBB = userInputs.get_model_constant_double("DBB");
    double DBV = userInputs.get_model_constant_double("DBV");
    double DVB = userInputs.get_model_constant_double("DVB");
    double DVV = userInputs.get_model_constant_double("DVV");
    double VM = userInputs.get_model_constant_double("VM");
    double zM = userInputs.get_model_constant_double("zM");
    double beta = userInputs.get_model_constant_double("beta");


    double delta = sqrt(2.0*epssq/W); //Equilibrium interface half-width
    double Mconst   = 2.0 * delta;  // psi mobility coeff
    double RT = R*T;

    //coefficient for eta in the Butler-Volmer exponential
    double eta_coeff = zM*(1.0-beta)*FarC/(RT);

    double cA_eq = 1.0-cB_eq;
    //temperature-dependent part of muB
    double muB_term1 = (-8856.94+157.48*T-26.908*T*std::log(T)+0.00189435*std::pow(T,2)-1.47721E-06*std::pow(T,3)+139250*(1.0/T));
    double muB_coeff2 = (4300 - 8.9*T);
    double muB_coeff3 = (27000 - 13.8*T);
    double muBFar = muB_term1 + muB_coeff2*cA_eq + muB_coeff3*cA_eq*(2*cB_eq - cA_eq) + RT*(std::log(cB_eq) + 1.0);

  // ================================================================
};
