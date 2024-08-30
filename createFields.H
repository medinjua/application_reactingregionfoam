#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<fluidReactionThermo> pThermo(fluidReactionThermo::New(mesh));
fluidReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureReference pressureReference(p, pimple.dict(), false);

mesh.setFluxRequired(p.name());


Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidReactionThermophysicalTransportModel>
thermophysicalTransport
(
    fluidReactionThermophysicalTransportModel::New
    (
        turbulence(),
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<combustionModel> reaction(combustionModel::New(thermo, turbulence()));


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

// ****************************************************************************
// Additions for porous catalyst simulations
// ****************************************************************************

const word inertSpecie(thermo.properties().lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo.properties())
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

// Gas constant units in J/(kmol * K)
const dimensionedScalar RuSI  = dimensionedScalar(
                                    dimensionSet(1,2,-2,-1,-1,0,0), 
                                    8314.463
                                );

Info<< "\nAdding porosity" << endl;
volScalarField epspore
(
    IOobject
    (
        "epspore",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Heat release rate [J/m3/s]
volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Field for second auxiliary porosity []
volScalarField epscatalyst
(
    IOobject
    (
        "epscatalyst",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epscatalyst",dimensionSet(0,0,0,0,0,0, 0),0)
);

// Check whether solid properties were given
if( !readScalar(thermo.properties().lookup("kappaSolid")) || !readScalar(thermo.properties().lookup("rhoSolid"))
    || !readScalar(thermo.properties().lookup("cvSolid"))  )
{
    FatalIOErrorIn(args.executable().c_str(), thermo.properties())
        << "Solid properties (kappaSolid, rhoSolid, cvSolid) not found! (thermophysicalProperties) "
        << exit(FatalIOError);
}

scalar kappa_solid = readScalar(thermo.properties().lookup("kappaSolid"));
scalar cv_solid = readScalar(thermo.properties().lookup("kappaSolid"));
scalar rho_solid = readScalar(thermo.properties().lookup("rhoSolid"));

Info<< "\nCreating field with solid thermal conductivity" << endl;
volScalarField kappas
(
    IOobject
    (
        "kappas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappas",dimensionSet(1,1,-3,-1,0,0,0),kappa_solid)
);

Info<< "\nCreating field with solid density" << endl;
volScalarField rhos
(
    IOobject
    (
        "rhos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhos",dimensionSet(1,-3,0,0,0,0,0),rho_solid)
);

Info<< "\nCreating field with solid specific heat capacity" << endl;
volScalarField cpvs
(
    IOobject
    (
        "cpvs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cpvs",dimensionSet(0,2,-2,-1,0,0,0),cv_solid)
);

Info<< "\nCreating field for effective thermal conductivity" << endl;
volScalarField kappaeffbycp
(
    IOobject
    (
        "kappaeffbycp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappaeffbycp",dimensionSet(1,-1,-1,0,0,0,0),0)
);