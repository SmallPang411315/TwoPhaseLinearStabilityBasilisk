double value_Heaviside(double phi, double Delta) {
	if(phi > Delta)
	    return 1.;
	else if(phi < -Delta)
		return 0.;
		else 
			return ( 1 + phi/Delta + sin( pi*phi/Delta )/pi )/2.;
}

double value_dHeaviside(double phi, double Delta) {
	if(phi > Delta)
	    return 0.;
	else if(phi < -Delta)
		return 0.;
		else 
			return cos(pi*phi/Delta)/2./Delta + 0.5/Delta;
}

double value_Heaviside_1(double phi, double Delta) {
			return (0.5+0.5*Delta*phi/sqrt(1+sq(Delta)*sq(phi)));
}

double value_dHeaviside_1(double phi, double Delta) {
			return (Delta/(2*sqrt(sq(Delta)*sq(phi) + 1)) - (cube(Delta)*sq(phi))/(2*sqrt(cube(sq(Delta)*sq(phi) + 1))));
}
