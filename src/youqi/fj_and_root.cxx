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
	UInt_t ntrials, evid, ncharged, nneutral, ncharged_g, nneutral_g;
	Float_t xsec, x, y, Q2, W2;
	Float_t e_jet, pt_jet, eta_jet, phi_jet, p_jet, theta_jet;
	Float_t e_jet_g, pt_jet_g, eta_jet_g, phi_jet_g, p_jet_g, theta_jet_g, delta_R, zg;
	Float_t e_quark, pt_quark, eta_quark, phi_quark, p_quark, theta_quark;
	Float_t e_electron, pt_electron, eta_electron, phi_electron, p_electron, theta_electron;
	Float_t e_photon, pt_photon, eta_photon, phi_photon;

	// initialize TTree
	TTree *tree = new TTree("Tree", "Tree");
	tree->Branch("ntrials", &ntrials, "ntrials/I");
	tree->Branch("evid", &evid, "evid/I");
	tree->Branch("ncharged", &ncharged, "ncharged/I");
	tree->Branch("nneutral", &nneutral, "nneutral/I");
	tree->Branch("ncharged_g", &ncharged_g, "ncharged_g/I");
	tree->Branch("nneutral_g", &nneutral_g, "nneutral_g/I");
	tree->Branch("xsec", &xsec, "xsec/F");
	tree->Branch("x", &x, "x/F");
	tree->Branch("y", &y, "y/F");
	tree->Branch("Q2", &Q2, "Q2/F");
	tree->Branch("W2", &W2, "W2/F");
	tree->Branch("e_jet", &e_jet, "e_jet/F");
	tree->Branch("pt_jet", &pt_jet, "pt_jet/F");
	tree->Branch("eta_jet", &eta_jet, "eta_jet/F");
	tree->Branch("phi_jet", &phi_jet, "phi_jet/F");
	tree->Branch("p_jet", &p_jet, "p_jet/F");
	tree->Branch("theta_jet", &theta_jet, "theta_jet/F");
	tree->Branch("e_jet_g", &e_jet_g, "e_jet_g/F");
	tree->Branch("pt_jet_g", &pt_jet_g, "pt_jet_g/F");
	tree->Branch("eta_jet_g", &eta_jet_g, "eta_jet_g/F");
	tree->Branch("phi_jet_g", &phi_jet_g, "phi_jet_g/F");
	tree->Branch("p_jet_g", &p_jet_g, "p_jet_g/F");
	tree->Branch("theta_jet_g", &theta_jet_g, "theta_jet_g/F");
	tree->Branch("delta_R", &delta_R, "delta_R/F");
	tree->Branch("zg",&zg,"zg/F");
	tree->Branch("e_quark", &e_quark, "e_quark/F");
	tree->Branch("pt_quark", &pt_quark, "pt_quark/F");
	tree->Branch("eta_quark", &eta_quark, "eta_quark/F");
	tree->Branch("phi_quark", &phi_quark, "phi_quark/F");
	tree->Branch("p_quark", &p_quark, "p_quark/F");
	tree->Branch("theta_quark", &theta_quark, "theta_quark/F");
	tree->Branch("e_electron", &e_electron, "e_electron/F");
	tree->Branch("pt_electron", &pt_electron, "pt_electron/F");
	tree->Branch("eta_electron", &eta_electron, "eta_electron/F");
	tree->Branch("phi_electron", &phi_electron, "phi_electron/F");
	tree->Branch("p_electron", &p_electron, "p_electron/F");
	tree->Branch("theta_electron", &theta_electron, "theta_electron/F");
	tree->Branch("e_photon", &e_photon, "e_photon/F");
	tree->Branch("pt_photon", &pt_photon, "pt_photon/F");
	tree->Branch("eta_photon", &eta_photon, "eta_photon/F");
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

		// Get struck quark index.
		int q;
		for (int i = 0; i < event.size(); i++)
		{
			if (event[i].status() == -23 && event[i].id() != 11)
			{
				q = i;
				break;
			}
		}

		// Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
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

		// run jet finding
		fj::JetDefinition jet_def(fj::antikt_algorithm, jetR, fj::WTA_pt_scheme);
		fj::ClusterSequence ca(parts_selected, jet_def);
		std::vector<fj::PseudoJet> jets_inclusive = ca.inclusive_jets();
		std::vector<fj::PseudoJet> jets = jetSelector(jets_inclusive);

		// soft drop jets
		std::vector<fj::PseudoJet> sdjets = FJUtils::soft_drop_jets(jets, 0.1, 0.0, jetR);

		// write jet properties to a TTree
		for (unsigned int ij = 0; ij < jets.size(); ij++)
		{
			// jet constituents loop
			for (unsigned int i = 0; i < jets[ij].constituents().size(); i++)
			{
				auto _part = fastjet::sorted_by_pt(jets[ij].constituents())[i];
				Pythia8::Particle *_part_py = _part.user_info<FJUtils::PythiaUserInfo>().getParticle();
				if (_part_py->isCharged() == true)
					ncharged += 1;
				else
					nneutral += 1;
			}
			// sdjet constituents loop
			for (unsigned int i = 0; i < sdjets[ij].constituents().size(); i++)
			{
				auto _part_g = fastjet::sorted_by_pt(sdjets[ij].constituents())[i];
				Pythia8::Particle *_part_py_g = _part_g.user_info<FJUtils::PythiaUserInfo>().getParticle();
				if (_part_py_g->isCharged() == true)
					ncharged_g += 1;
				else
					nneutral_g += 1;
			}

			evid = iEvent;
			xsec = pythia.info.sigmaGen();
			ntrials = pythia.info.nTried();
			e_jet = jets[ij].e();
			pt_jet = jets[ij].perp();
			eta_jet = jets[ij].eta();
			phi_jet = jets[ij].phi();
			p_jet = jets[ij].perp() * TMath::CosH(jets[ij].eta());
			theta_jet = jets[ij].theta();
			e_jet_g = sdjets[ij].e();
			pt_jet_g = sdjets[ij].perp();
			eta_jet_g = sdjets[ij].eta();
			phi_jet_g = sdjets[ij].phi();
			p_jet_g = sdjets[ij].perp() * TMath::CosH(sdjets[ij].eta());
			theta_jet_g = sdjets[ij].theta();
			delta_R = sdjets[ij].structure_of<fj::contrib::SoftDrop>().delta_R();
			zg = sdjets[ij].structure_of<fj::contrib::SoftDrop>().symmetry();
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
			pt_photon = pPhoton.pT();
			eta_photon = pPhoton.eta();
			phi_photon = pPhoton.phi();
			
			tree->Fill();

			ncharged = 0;
			nneutral = 0;
			ncharged_g = 0;
			nneutral_g = 0;
	
		} //jet
	}	 //event

	// write and close the output file
	fout.Write();
	fout.Close();
	cout << "[i] file written: " << fout.GetName() << endl;
	// Done.
	return 0;
} // fj_and_root
} // namespace youqi
