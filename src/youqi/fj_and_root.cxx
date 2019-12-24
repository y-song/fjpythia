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
#include <math.h>

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
	Float_t phi_photon;

	// new particle variables
	std::vector<float> pt_jet;
	std::vector<float> phi_jet;
	std::vector<float> z;
	std::vector<float> j;

	// initialize TTree
	TTree *tree1 = new TTree("Tree1", "Tree1");
	tree1->Branch("ntrials", &ntrials, "ntrials/I");
	tree1->Branch("evid", &evid, "evid/I");
	tree1->Branch("xsec", &xsec, "xsec/F");
	tree1->Branch("x", &x, "x/F");
	tree1->Branch("y", &y, "y/F");
	tree1->Branch("Q2", &Q2, "Q2/F");
	tree1->Branch("W2", &W2, "W2/F");
	tree1->Branch("phi_photon", &phi_photon, "phi_photon/F");
	tree1->Branch("pt_jet", &pt_jet);
	tree1->Branch("phi_jet", &phi_jet);
	tree1->Branch("z", &z);
	tree1->Branch("j", &j);

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
		
	        z.clear();
       		j.clear();
		pt_jet.clear();
		phi_jet.clear();

		// four-momenta of proton, electron, virtual photon/Z^0/W^+-.
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

		evid = iEvent;
		xsec = pythia.info.sigmaGen();
		ntrials = pythia.info.nTried();
		phi_photon = pPhoton.phi();

		// run jet finding
		fj::JetDefinition jet_def(fj::antikt_algorithm, jetR);
		fj::ClusterSequence ca(parts_selected, jet_def);
		std::vector<fj::PseudoJet> jets_inclusive = ca.inclusive_jets();
		std::vector<fj::PseudoJet> jets = jetSelector(jets_inclusive);

		// loop over jets
		for (unsigned int ij = 0; ij < jets.size(); ij++)
		{
			// loop over particles
			for (int i = 0; i < jets[ij].constituents().size(); i++)
			{
				Pythia8::Particle *_p = jets[ij].constituents()[i].user_info<FJUtils::PythiaUserInfo>().getParticle();
               			// save only charged hadrons
                		if (_p->isCharged() && _p->isHadron())
                		{
               	     			double pxj, pyj, pzj, pxh, pyh, pzh, p_jet, cross;
					pxj = jets[ij].px();
                    			pyj = jets[ij].py();
                    			pzj = jets[ij].pz();
                    			pxh = _p->px();
                    			pyh = _p->py();
                    			pzh = _p->pz();
                    			p_jet = sqrt(pxj*pxj + pyj*pyj + pzj*pzj);
                   			cross = sqrt( pow((pyj*pzh-pyh*pzj),2.0) + pow((pxj*pzh-pzj*pxh),2.0) + pow((pxj*pyh-pyj*pxh),2.0) );
                    			z.push_back((pxj*pxh + pyj*pyh + pzj*pzh) / (p_jet*p_jet));
                    			j.push_back(cross / p_jet);
				    	pt_jet.push_back(jets[ij].perp());
        			    	phi_jet.push_back(jets[ij].phi());
                		}
			}	
		}	
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
