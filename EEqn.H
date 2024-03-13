{
    
    if(!isothermalFluid){
        //**********************************************************************
        // Modification catalyst simulations
        // TODO: Modify alphaeff.H here
        #include "alphaeff.H"
        // Reinitialize Qdot to zero
        Qdot *= 0.0;

        // First, calculate the heat release rate by the chemistry pointer
        // Note the HRR units (J/m3/s)
        // Note: a call to chemistry.calculate() must have been done before
        volScalarField::Internal HRR(chemistry.Qdot().ref());
        // Now set the HRR in the catalyst cells
        forAll(cellsCatalyst, cc)
        {
            label celli = cellsCatalyst[cc];
            Qdot[celli] = HRR[celli]*(1.0 - epspore[celli]);
        }
        //******************************************************************

        volScalarField& he = thermo.he();

        fvScalarMatrix EEqn
        (
            fvm::ddt((epspore*rho*cp_f + (1.0-epspore)*rhos*cpvs)/cp_f, he) + epspore*mvConvection->fvmDiv(phi, he)
          //+ epspore*fvc::ddt(rho, K) + fvc::div(phi, K)
          + (
                he.name() == "e"
                ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
                : -dpdt
            )
          //+ thermophysicalTransport->divq(he) // replaced to allow kappaeff
                                                // ->implicit: fvm::laplacian
                                                // ->explicit: fvc::laplacian
          - fvm::laplacian(kappaeffbycp, he)    // see definition of
                                                // thermophysicalTransport>divq
          ==
            //reaction->Qdot()          // commented out for replacement *******
            Qdot                        // replaced HRR             ************
          + fvModels.source(rho, he)
        );

        EEqn.relax();

        fvConstraints.constrain(EEqn);

        EEqn.solve();

        fvConstraints.constrain(he);

        thermo.correct();

    }

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
