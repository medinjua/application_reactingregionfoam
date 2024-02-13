// Overwrite alphaeff tensor to the corresponding 
// fluid thermal diffusivity for zones without catalyst
volScalarField kappa_f = thermo.kappa();
volScalarField cp_f    = thermo.Cp();

forAll(cellsNoCatalyst, cc)
{
    label celli          = cellsNoCatalyst[cc];
    kappaeffbycp[celli]  = kappa_f[celli] / cp_f[celli];
}

// The effective thermal conductivity is calculated by 
// a weighted average of solid and fluid thermal conductivities
// That is, keff = ksolid*(1-epsilon) + kfluid*epsilon
forAll(cellsCatalyst, cc)
{
    label celli          = cellsCatalyst[cc];
    kappaeffbycp[celli]  = (kappa_f[celli]*epspore[celli] + (1.0 - epspore[celli])*kappas[celli])/ cp_f[celli];
}