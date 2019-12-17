//- relative permeability (kr)
krModel->correct();
const volScalarField& kr = krModel->kr();
const volScalarField& dkrdS = krModel->dkrdS();

