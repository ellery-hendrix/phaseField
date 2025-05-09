// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for
// each function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void
variableAttributeLoader::loadVariableAttributes()
{

    // Variable 0: domain parameter for grain
    set_variable_name               (0,"psi");
    set_variable_type               (0,SCALAR);
    set_variable_equation_type      (0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "psi");
    set_dependencies_gradient_term_RHS(0, "psi, grad(mupsi)");

    // Variable 1: chemical potential for domain parameter
    set_variable_name               (1,"mupsi");
    set_variable_type               (1,SCALAR);
    set_variable_equation_type      (1,AUXILIARY);

    set_dependencies_value_term_RHS(1, "psi");
    set_dependencies_gradient_term_RHS(1, "grad(psi)");

    // Variable 2: concentration of B, fast diffuser which also leaches out
    set_variable_name               (2,"cB");
    set_variable_type               (2,SCALAR);
    set_variable_equation_type      (2,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(2, "cB, cV, psi, grad(psi)");
    set_dependencies_gradient_term_RHS(2, "grad(cB),grad(cV)");

    // Variable 3: concentration of vacancies
    set_variable_name               (3,"cV");
    set_variable_type               (3,SCALAR);
    set_variable_equation_type      (3,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(3, "cB, cV, psi, grad(psi)");
    set_dependencies_gradient_term_RHS(3, "grad(cV),grad(cB)");
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time
// dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a
// list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one
// proportional to the test function and one proportional to the gradient of the
// test function. The index for each variable in this list corresponds to the
// index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{

// --- Getting the values and derivatives of the model variables ---

// Timestep
scalarvalueType delt = constV(userInputs.dtValue);
// psi
scalarvalueType psi = variable_list.get_scalar_value(0);
scalargradType gradpsi = variable_list.get_scalar_gradient(0);

// mupsi
scalargradType gradmupsi = variable_list.get_scalar_gradient(1);

// cB
scalarvalueType cB = variable_list.get_scalar_value(2);
scalargradType gradcB = variable_list.get_scalar_gradient(2);

// cV
scalarvalueType cV = variable_list.get_scalar_value(3);
scalargradType gradcV = variable_list.get_scalar_gradient(3);

// --- Calculation of terms needed in multiple expressions  ---
scalarvalueType invpsicp = constV(0.0);
scalarvalueType invpsi = constV(0.0);

for (unsigned int j=0; j<psi.size();j++){
    invpsicp[j]=psi[j];
    if (psi[j] < lthresh)
        invpsicp[j] = 1.0/lthresh;
    if (psi[j] >= lthresh)
        invpsicp[j] = 1.0/psi[j];
}
invpsi = invpsicp;

// --- Setting the expressions for the terms in the governing equations ---
scalarvalueType Mpsi = Mconst * psi;


// Magnitude of the gradient of the domain parameter
scalarvalueType maggradpsi = constV(0.0);
// Dot products
scalarvalueType gradpsixcB = constV(0.0);
scalarvalueType gradpsixcV = constV(0.0);
for (int k=0; k<dim; k++){
    // magnitude of gradient psi
    maggradpsi = maggradpsi+ gradpsi[k]*gradpsi[k];
    // dot product between grad(psi) and grad(cB)
    gradpsixcB = gradpsixcB + gradpsi[k]*gradcB[k];
    // dot product between grad(psi) and grad(cB)
    gradpsixcV = gradpsixcV + gradpsi[k]*gradcV[k];

}
maggradpsi = std::sqrt(maggradpsi);

// chemical potential of B
scalarvalueType cA = constV(1.0)-cB;
scalarvalueType muB = constV(muB_term1) + constV(muB_coeff2)*cA + constV(muB_coeff3)*cA*(constV(2.0)*cB - cA) + constV(RT)*(std::log(cB) + constV(1.0));

//Reaction Current
scalarvalueType irxni = constV(0.0);    //reaction current for a single grain i
scalarvalueType irxn = constV(0.0);     //reaction current for the whole domain
scalarvalueType eta=(muB-constV(muBFar))/constV(FarC);
irxni=icorr*std::exp(constV(eta_coeff)*eta);

// --- Residuals ---
scalarvalueType rpsi  = psi;
scalargradType rgradpsi = -Mpsi * delt * gradmupsi;

scalarvalueType rcB = cB - delt*constV(DBB)*invpsi*gradpsixcB - delt*constV(DBV)*invpsi*gradpsixcV - delt*irxn*constV(VM/(zM*FarC))*invpsi*maggradpsi;
scalargradType rgradcB = -delt*(DBB*gradcB)-delt*(DBV*gradcV);

scalarvalueType rcV = cV - delt*constV(DVB)*invpsi*gradpsixcB - delt*constV(DVV)*invpsi*gradpsixcV;
scalargradType rgradcV = -delt*(DVB*gradcV)-delt*(DVV*gradcV);

// --- Submitting the terms for the governing equations ---

variable_list.set_scalar_value_term_RHS(0,rpsi);
variable_list.set_scalar_gradient_term_RHS(0,rgradpsi);

variable_list.set_scalar_value_term_RHS(2,rcB);
variable_list.set_scalar_gradient_term_RHS(2,rgradcB);

variable_list.set_scalar_value_term_RHS(3,rcV);
variable_list.set_scalar_gradient_term_RHS(3,rgradcV);


}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time
// independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are
// not explicit time-dependent equations. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is
// given by "q_point_loc". The function outputs two terms to variable_list --
// one proportional to the test function and one proportional to the gradient of
// the test function. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. The
// (x,y,z) location of that quadrature point is given by "q_point_loc". The
// function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function -- for the
// left-hand-side of the equation. The index for each variable in this list
// corresponds to the index given at the top of this file. If there are multiple
// elliptic equations, conditional statements should be sed to ensure that the
// correct residual is being submitted. The index of the field being solved can
// be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{}
