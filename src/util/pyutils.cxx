#include <fjpythia/util/pyutils.h>
#include <fjpythia/util/argparser.h>
#include <fjpythia/util/strutil.h>

#include <Pythia8/Pythia.h>
#include <TString.h>

#include <iostream>
using namespace std;

namespace PythiaUtils
{
	std::vector<int> find_outgoing_hard_electrons(Pythia8::Pythia *pythia)
	{
		std::vector<int> v;
		for (unsigned int i = 0; i < pythia->event.size(); i++)
		{
			if (!pythia->event[i].isFinal()) continue;
			if (pythia->event[i].id() == 11)
			{
				if (pythia->event[i].status() == 23)
					v.push_back(i);
			}
		}
		return v;
	}

	void cook_pythia_settings(Pythia8::Pythia *pythia)
	{
		auto &args = FJPyUtil::ArgParser::Instance();

		// settings from a cmnd file?
		if (args.isSet("--pythia-config"))
		{
			pythia->readFile(args.getOpt("--pythia-config").c_str());
			args.addOpts("--pythia-process-configured");
		}

		// settings from a command line?
		std::string pythiaOpt = args.getOpt("--pythia");
		auto pyopts = StrUtil::split_to_vector(pythiaOpt.c_str(), ",");
		for (auto o : pyopts)
		{
			StrUtil::replace_substring(o, "_", " ");
			pythia->readString(o.c_str());
			args.addOpts("--pythia-process-configured");
		}

		int nDebugEvents = args.getOptInt("--debug-events", 0);
		pythia->readString(TString::Format("Next:numberShowEvent=%d", nDebugEvents).Data());
		pythia->readString(TString::Format("Next:numberShowInfo=%d", nDebugEvents).Data());
		pythia->readString(TString::Format("Next:numberShowProcess=%d", nDebugEvents).Data());
		pythia->readString(TString::Format("Next:numberCount=%d", nDebugEvents).Data());

		if (args.isSet("--hardQCD"))
		{
			pythia->readString("HardQCD:all=on");
			args.addOpts("--pythia-process-configured");
		}

		if (args.isSet("--eic-dis") or args.isSet("--eic-lowQ2"))
		{
			pythia->readString("Beams:idA=11");
			pythia->readString("Beams:idB=2212");
			pythia->readString("Beams:eA=20");
			pythia->readString("Beams:eB=250");
			pythia->readString("Beams:frameType=2");
			pythia->readString("Init:showChangedSettings=on");
			pythia->readString("Main:timesAllowErrors=10000");
			if (args.isSet("--eic-dis"))
			{
				pythia->readString("WeakBosonExchange:ff2ff(t:gmZ)=on");
				pythia->readString("PhaseSpace:Q2Min=10");
				pythia->readString("SpaceShower:pTmaxMatch=2");
				pythia->readString("PDF:lepton=off");
				pythia->readString("TimeShower:QEDshowerByL=off");

				args.addOpts("--pythia-process-configured");
			}
			if (args.isSet("--eic-lowQ2"))
			{
				pythia->readString("HardQCD:all=on");
				pythia->readString("PDF:lepton2gamma=on");
				pythia->readString("Photon:Q2max=1.");
				pythia->readString("Photon:Wmin=10.");
				pythia->readString("PhaseSpace:pTHatMin=2.");
				pythia->readString("PhotonParton:all=on");
				pythia->readString("Photon:ProcessType=0");

				args.addOpts("--pythia-process-configured");
			}
		}

		if (!args.isSet("--pythia-process-configured"))
		{
			cout << "[w] this will likely fail. enable a process with a flag or --pythia list,of,settings,..." << endl;
			args.addOpts("--debug");
		}

	}
}