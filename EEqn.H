{
    
    if(!isothermalFluid){
        //**********************************************************************
        // Modification catalyst simulations
        #include "alphaeff.H"
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
            (1.0 - epspore)*reaction->Qdot()
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
