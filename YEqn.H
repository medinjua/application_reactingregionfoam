tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

reaction->correct();

forAll(Y, i)
{
    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(epspore*rho, Yi)
          + epspore*mvConvection->fvmDiv(phi, Yi)
          // Note: default laminar transport model for species diffusion 
          // is unityLewisFourier, the species diffusion coefficient is 
          // then simply the thermal diffusion coefficient
          // note that the documentation for divj says
          // -fvm::laplacian(this->alpha()*this->DEff(Yi), Yi)
          // this->alpha() refers to the alphaField, simply a volume fraction
          // DEff is simply the thermal diffusion coefficient
          + epspore*thermophysicalTransport->divj(Yi)
         ==
            (1.0-epspore)*reaction->R(Yi)
          + fvModels.source(rho, Yi)
        );

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);
    }
}

composition.normalise();
