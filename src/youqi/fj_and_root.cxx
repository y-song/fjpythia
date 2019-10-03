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
	UInt_t ntrials, evid, id_quark;
	Float_t xsec, x, y, Q2, W2;
	Float_t e_quark, pt_quark, eta_quark, phi_quark, p_quark, theta_quark;
	Float_t e_electron, pt_electron, eta_electron, phi_electron, p_electron, theta_electron;
	Float_t e_photon, eta_photon, phi_photon;

	// new particle variables
	std::vector<int> jetid;
	std::vector<float> pt_jet;
	std::vector<float> eta_jet;
	std::vector<int> id_lead_part;
	std::vector<int> id_part;
	std::vector<float> pT_part;
	std::vector<float> p_part;
	std::vector<float> eta_part;
	std::vector<float> theta_part;
	std::vector<float> y_part;
	std::vector<float> phi_part;
	std::vector<float> e_part;

	// initialize TTree
	TTree *tree1 = new TTree("Tree1", "Tree1");
	tree1->Branch("ntrials", &ntrials, "ntrials/I");
	tree1->Branch("evid", &evid, "evid/I");
	tree1->Branch("xsec", &xsec, "xsec/F");
	tree1->Branch("x", &x, "x/F");
	tree1->Branch("y", &y, "y/F");
	tree1->Branch("Q2", &Q2, "Q2/F");
	tree1->Branch("W2", &W2, "W2/F");
	tree1->Branch("id_quark", &id_quark, "id_quark/I");
	tree1->Branch("e_quark", &e_quark, "e_quark/F");
	tree1->Branch("pt_quark", &pt_quark, "pt_quark/F");
	tree1->Branch("eta_quark", &eta_quark, "eta_quark/F");
	tree1->Branch("phi_quark", &phi_quark, "phi_quark/F");
	tree1->Branch("p_quark", &p_quark, "p_quark/F");
	tree1->Branch("theta_quark", &theta_quark, "theta_quark/F");
	tree1->Branch("e_electron", &e_electron, "e_electron/F");
	tree1->Branch("pt_electron", &pt_electron, "pt_electron/F");
	tree1->Branch("eta_electron", &eta_electron, "eta_electron/F");
	tree1->Branch("phi_electron", &phi_electron, "phi_electron/F");
	tree1->Branch("p_electron", &p_electron, "p_electron/F");
	tree1->Branch("theta_electron", &theta_electron, "theta_electron/F");
	tree1->Branch("e_photon", &e_photon, "e_photon/F");
	tree1->Branch("jetid", &jetid);
	tree1->Branch("pt_jet", &pt_jet);
	tree1->Branch("eta_jet", &eta_jet);
	tree1->Branch("id_lead_part", &id_lead_part);
	tree1->Branch("id_part", &id_part);
	tree1->Branch("pT_part", &pT_part);
	tree1->Branch("p_part", &p_part);
	tree1->Branch("eta_part", &eta_part);
	tree1->Branch("theta_part", &theta_part);
	tree1->Branch("phi_part", &phi_part);
	tree1->Branch("e_part", &e_part);

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
	int nEv = args.getOptInt("--nev", 4); // default will be 4 event(!)
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

	// Begin event loop. Generate event. Skip if error..
	int num_jet = 0;
	LoopUtil::TPbar pbar(nEv);
	for (int iEvent = 0; iEvent < nEv; ++iEvent)
	{
		pbar.Update();
		if (!pythia.next())
			continue;
		
		jetid.clear();
		pt_jet.clear();
		eta_jet.clear();
		id_lead_part.clear();
		id_part.clear();
		pT_part.clear();
		p_part.clear();
		eta_part.clear();
		theta_part.clear();
		phi_part.clear();
		e_part.clear();

		// get struck quark index
		int q;
		for (int i = 0; i < event.size(); i++)
		{
			if (event[i].status() == -23 && event[i].id() != 11)
			{
				q = i;
				break;
			}
		}

		// four-momenta of proton, electron, virtual photon/Z^0/W^+-.
		Vec4 pProton = event[1].p();
		Vec4 peIn = event[4].p();
		Vec4 peOut = event[6].p();
		Vec4 pPhoton = peIn - peOut;
		Vec4 pQuark = event[q].p();

		// Q2, W2, Bjorken x, y, nu.
		Q2 = -pPhoton.m2Calc();
		W2 = (pProton + pPhoton).m2Calc();
		x = Q2 / (2. * pProton * pPhoton);
		y = (pProton * pPhoton) / (pProton * peIn);

		auto parts = FJUtils::getPseudoJetsFromPythia(&pythia, true); // only_final==true
		std::vector<fj::PseudoJet> parts_selected = partSelector(parts);

		evid = iEvent;
		xsec = pythia.info.sigmaGen();
		ntrials = pythia.info.nTried();
		id_quark = event[q].id();
		e_quark = event[q].e();
		pt_quark = event[q].pT();
		eta_quark = event[q].eta();
		phi_quark = event[q].phi();
		p_quark = event[q].pT() * TMath::CosH(event[q].eta());
		theta_quark = event[q].theta();
		e_electron = event[6].e();
		pt_electron = event[6].pT();
		eta_electron = event[6].eta();
		phi_electron = event[6].phi();
		p_electron = event[6].pT() * TMath::CosH(event[6].eta());
		theta_electron = event[6].theta();
		e_photon = pPhoton[0];

		// run jet finding
		fj::JetDefinition jet_def(fj::antikt_algorithm, jetR);
		fj::ClusterSequence ca(parts_selected, jet_def);
		std::vector<fj::PseudoJet> jets_inclusive = ca.inclusive_jets();
		std::vector<fj::PseudoJet> jets = jetSelector(jets_inclusive);

		// save jet kinematics: loop over jets
		for (unsigned int ij = 0; ij < jets.size(); ij++)
		{
			// loop over particles
			for (int i = 0; i < jets[ij].constituents().size(); i++)
			{
				jetid.push_back(num_jet + ij);
				pt_jet.push_back(jets[ij].perp());
				eta_jet.push_back(jets[ij].eta());

				auto _lead = fastjet::sorted_by_pt(jets[ij].constituents())[0];
				Pythia8::Particle *_lead_p = _lead.user_info<FJUtils::PythiaUserInfo>().getParticle();
				id_lead_part.push_back(_lead_p->id());

				Pythia8::Particle *_p = jets[ij].constituents()[i].user_info<FJUtils::PythiaUserInfo>().getParticle();
				id_part.push_back(_p->id());
				pT_part.push_back(jets[ij].constituents()[i].perp());
				p_part.push_back(jets[ij].constituents()[i].perp() * TMath::CosH(jets[ij].constituents()[i].eta()));
				eta_part.push_back(jets[ij].constituents()[i].eta());
				theta_part.push_back(jets[ij].constituents()[i].theta());
				phi_part.push_back(jets[ij].constituents()[i].phi());
				e_part.push_back(jets[ij].constituents()[i].e());
			}
		}
		num_jet = num_jet + jets.size();
		tree1->Fill();

	} //event

	// write and close the output file
	fout.Write();
	fout.Close();
	cout << "[i] file written: " << fout.GetName() << endl;
	// Done.
	return 0;
} // fj_and_root
} // namespace youqi
