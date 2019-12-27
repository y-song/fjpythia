#include <fjpythia/eic/example.h>
#include <fjpythia/util/argparser.h>
#include <fjpythia/util/strutil.h>
#include <fjpythia/util/fjutils.h>
#include <fjpythia/util/looputil.h>
#include <fjpythia/util/pyutils.h>

#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>

#include <fastjet/Selector.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/contrib/SoftDrop.hh>

namespace fj = fastjet;

#include "Pythia8/Pythia.h"
using namespace Pythia8;

namespace youqi
{
int fj_and_root()
{
	auto &args = FJPyUtil::ArgParser::Instance();
	args.dump();

	// open an output file
	std::string foutname = args.getOpt("--out", "default_output.root");
	cout << "[i] output file name: " << foutname << endl;
	TFile fout(foutname.c_str(), "recreate");
	fout.cd();

	// new branches
	UInt_t ntrials, evid;
	Float_t xsec, x, y, Q2, W2;
	Float_t pt_jet, phi_jet;
	Float_t phi_photon;

	// initialize TTree
	TTree *tree = new TTree("Tree", "Tree");
	tree->Branch("ntrials", &ntrials, "ntrials/I");
	tree->Branch("evid", &evid, "evid/I");
	tree->Branch("xsec", &xsec, "xsec/F");
	tree->Branch("x", &x, "x/F");
	tree->Branch("y", &y, "y/F");
	tree->Branch("Q2", &Q2, "Q2/F");
	tree->Branch("W2", &W2, "W2/F");
	tree->Branch("pt_jet", &pt_jet, "pt_jet/F");
	tree->Branch("phi_jet", &phi_jet, "phi_jet/F");
	tree->Branch("phi_photon", &phi_photon, "phi_photon/F");

	// intialize PYTHIA
	Pythia pythia;
	Event &event = pythia.event;
	PythiaUtils::cook_pythia_settings(&pythia);
	if (!pythia.init())
	{
		cout << "[e] pythia init failed." << endl;
		return -1;
	}

	// generate and analyze events
	int nEv = args.getOptInt("--nev", 1); // default will be 1 event(!)
	double jetR = args.getOptDouble("--jetR", 1.0);
	double minJetPt = args.getOptDouble("--minJetPt", 0.0);
	double maxPartEta = std::abs(args.getOptDouble("--maxParticleEta", 4.5));
	double minPartPt = args.getOptDouble("--minPartPt", 0.25);

	// dump some parameters of the analysis
	cout << "[i] configuration: " << endl
		 << "    events:        " << nEv << endl
		 << "    jetR:          " << jetR << endl
		 << "    minJetPt:      " << minJetPt << endl
		 << "    maxPartEta:    " << maxPartEta << endl
		 << "    output:        " << foutname << endl;

	fj::Selector partSelector = fastjet::SelectorAbsEtaMax(maxPartEta) * fastjet::SelectorPtMin(minPartPt);
	fj::Selector jetSelector = fastjet::SelectorAbsEtaMax(maxPartEta - jetR - 0.01) * fastjet::SelectorPtMin(minJetPt);

	// Begin event loop. Generate event. Skip if error. List first few.
	LoopUtil::TPbar pbar(nEv);
	for (int iEvent = 0; iEvent < nEv; ++iEvent)
	{
		pbar.Update();
		if (!pythia.next())
			continue;

		// Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
		Vec4 pProton = event[1].p();
		Vec4 peIn = event[4].p();
		Vec4 peOut = event[6].p();
		Vec4 pPhoton = peIn - peOut;
		
		// Q2, W2, Bjorken x, y, nu.
		Q2 = -pPhoton.m2Calc();
		W2 = (pProton + pPhoton).m2Calc();
		x = Q2 / (2. * pProton * pPhoton);
		y = (pProton * pPhoton) / (pProton * peIn);

		auto parts = FJUtils::getPseudoJetsFromPythia(&pythia, true); // only_final==true
		std::vector<fj::PseudoJet> parts_selected = partSelector(parts);

		// run jet finding
		fj::JetDefinition jet_def(fj::antikt_algorithm, jetR);
		fj::ClusterSequence ca(parts_selected, jet_def);
		std::vector<fj::PseudoJet> jets_inclusive = ca.inclusive_jets();
		std::vector<fj::PseudoJet> jets = jetSelector(jets_inclusive);

		// write jet properties to a TTree
		for (unsigned int ij = 0; ij < jets.size(); ij++)
		{

			evid = iEvent;
			xsec = pythia.info.sigmaGen();
			ntrials = pythia.info.nTried();
			pt_jet = jets[ij].perp();
			phi_jet = jets[ij].phi();
			phi_photon = pPhoton.phi();
			
			tree->Fill();
	
		} //jet
		
	} //event

	// write and close the output file
	fout.Write();
	fout.Close();
	cout << "[i] file written: " << fout.GetName() << endl;
	// Done.
	return 0;
} // fj_and_root
} // namespace youqi
