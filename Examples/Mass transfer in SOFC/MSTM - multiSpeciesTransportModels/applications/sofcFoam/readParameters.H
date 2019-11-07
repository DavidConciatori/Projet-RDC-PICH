    Info<< "Reading electrolyteDictionary\n" << endl;
    
    IOdictionary electrolyteDictionary
    (
        IOobject
        (
            "electrolyte",
            regions.mesh("electrolyte").time().constant(),
            regions.mesh("electrolyte"),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    
    Info<< "Reading solutionSwitch\n" << endl;

    IOdictionary solutionSwitch
    (
        IOobject
        (
            "solutionSwitch",
            "modelInput",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    
    word solverType(solutionSwitch.subDict("switch").lookup("solverType"));
    scalar startValue = 0;
    scalarList values;
    scalar perCentTopLevelTol = 0.1;
    scalar Eocv = readScalar(electrolyteDictionary.lookup("Eocv"));
    scalarList topLevelTollerance;
    label solverSwitch = 1;

    switch(solverType[0])
    {
        case 'i':
        case 'I':
        case 'c':
        case 'C':
            Info << "\n\nSolving in term of current density\n\n";
            startValue = electrolyteDictionary.lookupOrDefault<scalar>("startI", 0);
            values = scalarList(electrolyteDictionary.lookup("currentDensity"));
            topLevelTollerance = 1e-3 * values / values;
            solverSwitch = 1;
            break;
            
        case 'v':
        case 'V':
        case 'p':
        case 'P':
            Info << "\n\nSolving in term of voltage\n\n";
            startValue = electrolyteDictionary.lookupOrDefault<scalar>("startV", 0);
            values = scalarList(electrolyteDictionary.lookup("voltage"));
            topLevelTollerance = perCentTopLevelTol * (Eocv - values) / 5e-3;
            Info << "\n\n" << topLevelTollerance << "\n\n";
            solverSwitch = 0;
            break;
            
        default:
            Info << "\nInvalid solver type... Please select in solutionSwitch dictionary:"
                << "\n i, I, currentDensity or CurrentDensity to impose the current density value"
                << "\n v or V to impose the voltage value" << endl; 
            Foam::FatalError.exit();
            break;
    }
  
  
    const dictionary& switchDict = solutionSwitch.subDict("switch");

    bool massTransportSwitch = readScalar(switchDict.lookup("massTransport"));
    bool electricalSwitch = readScalar(switchDict.lookup("electrical"));
    bool printFlowSwitch = readScalar(switchDict.lookup("printFlow"));

    bool anodeSwitch = readScalar(switchDict.lookup("anode"));
    bool anodeChannelSwitch = readScalar(switchDict.lookup("anodeChannel"));
    bool cathodeSwitch = readScalar(switchDict.lookup("cathode"));
    bool cathodeChannelSwitch = readScalar(switchDict.lookup("cathodeChannel"));

    bool reactiveAnode = readScalar(switchDict.lookup("heterogeneousReactions"));
