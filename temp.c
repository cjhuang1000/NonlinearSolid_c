/* =========Integrand evaluation functions====== */
/* % Subcell elemental matrices are calculating from integral. Depending
% of the elementary functions involved (N_x^xi, N_y^xi or N^p), the
% formula is diffrent */
/*
double boundFunct_MS_alpha(unsigned int k, unsigned int l, double x, double y,double alpha, double dx){
	
	double a = dx^2 * boudfunc.coeff[1][k] * boudfunc.coeff[1][l]*( 
    (2*(x^3)/6 + (boudfunc.coeff(2,k)+boudfunc.coeff(2,l))*(x^2)/2 + boundFunct.coeff(2,k)*boundFunct.coeff(2,l)*x)   *   (2*(y^3)/6  + (boundFunct.coeff(3,k)+boundFunct.coeff(3,l))*(y^2)/2  + boundFunct.coeff(3,k)*boundFunct.coeff(3,l)* y      ) ...
                           + (-alpha)  *(2*(x^2)/2 + (boundFunct.coeff(2,k)+boundFunct.coeff(2,l))* x      + boundFunct.coeff(2,k)*boundFunct.coeff(2,l)  )   *   (2*(y^4)/24 + (boundFunct.coeff(3,k)+boundFunct.coeff(3,l))*(y^3)/6  + boundFunct.coeff(3,k)*boundFunct.coeff(3,l)*(y^2)/2 ) ...
                           + (-alpha)^2*(2* x      + (boundFunct.coeff(2,k)+boundFunct.coeff(2,l))                                                        )   *   (2*(y^5)/120+ (boundFunct.coeff(3,k)+boundFunct.coeff(3,l))*(y^4)/24 + boundFunct.coeff(3,k)*boundFunct.coeff(3,l)*(y^3)/6 ) ...
                           + (-alpha)^3*(2                                                                                                                )   *   (2*(y^6)/720+ (boundFunct.coeff(3,k)+boundFunct.coeff(3,l))*(y^5)/120+ boundFunct.coeff(3,k)*boundFunct.coeff(3,l)*(y^4)/24) ...
                        );
	return a;
}
*/