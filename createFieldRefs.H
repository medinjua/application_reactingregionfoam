const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();

// ****************************************************************************
// Additions for porous catalyst simulations
// ****************************************************************************

// chemistry model required
basicChemistryModel& chemistry = pChemistry();

// Define cell set list with and without catalyst zone
const word cellSetCatalyst = "CatalystRegionZone";
label zoneID = mesh.cellZones().findZoneID(cellSetCatalyst);
if (zoneID == -1)
{
    FatalErrorIn("createFieldRefs")
        << "Cannot find cellZone " << cellSetCatalyst << endl
        << "Valid cellZones are " << mesh.cellZones().names()
        << exit(FatalError);
}
const labelList& cellsCatalyst = mesh.cellZones()[zoneID];

const word cellSetNoCatalyst = "gas";
zoneID = mesh.cellZones().findZoneID(cellSetNoCatalyst);
if (zoneID == -1)
{
    FatalErrorIn("createFieldRefs")
        << "Cannot find cellZone " << cellSetNoCatalyst << endl
        << "Valid cellZones are " << mesh.cellZones().names()
        << exit(FatalError);
}
const labelList& cellsNoCatalyst = mesh.cellZones()[zoneID];

// Overwrite porosity (epsilon = 1) for zones without catalyst
forAll( cellsNoCatalyst, cc ){
    label celli    = cellsNoCatalyst[cc];
    epspore[celli] = 1.0;
}

// Check for option to run isothermal simulation (deactivates energy eqn)
bool isothermalFluid = readBool(thermo.properties().lookup("isothermalFlow"));
if( isothermalFluid ){
    Info<< "\nRunning isothermalFluid mode (no energy eqn solved)" << endl;
}