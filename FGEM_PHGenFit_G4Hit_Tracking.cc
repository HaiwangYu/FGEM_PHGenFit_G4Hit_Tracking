/*!
 *  \file		testPHGenFit.cc
 *  \brief		Program to demonstrate the usage of PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

//STL
#include <vector>

//BOOST
#include<boost/make_shared.hpp>

#define SMART(expr) boost::shared_ptr<expr>
#define NEW(expr) boost::make_shared<expr>

//sPHENIX
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>

//ROOT
#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TRandom.h>
#include <TClonesArray.h>
//#include <TColor.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

//GenFit
#include <GenFit/AbsTrackRep.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>

//PHGenFit
#include <phgenfit/Fitter.h>
#include <phgenfit/Track.h>
#include <phgenfit/Measurement.h>
#include <phgenfit/PlanarMeasurement.h>

#define LogPrint(exp)    std::cout<<__LINE__<<" : "<< #exp <<" : "<< exp <<"\n"

#define LogDEBUG(exp)    std::cout<<__LINE__<<" : DEBUG: "<< exp <<"\n"

#define WILD_NUM  -9999.;

double _FGEM_r_resolution = 1.; //1cm
double _FGEM_phi_resolution = 100E-4; //50um

//void pause() {
//  std::cout << "Press ENTER to continue..." << std::flush;
//  std::cin.clear();  // use only if a human is involved
//  std::cin.flush();  // use only if a human is involved
//  std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
//}

void get_seed(TVector3& seed_pos, TVector3& seed_mom, TMatrixDSym& seed_cov,
		PHG4Particle* particle = NULL, bool do_smearing = true) {

	seed_pos.SetXYZ(0, 0, 0);
	seed_mom.SetXYZ(0, 0, 10);
	seed_cov.ResizeTo(6, 6);

	for (int i = 0; i < 3; i++) {
		seed_cov[i][i] = _FGEM_phi_resolution * _FGEM_phi_resolution;
	}

	for (int i = 3; i < 6; i++) {
		seed_cov[i][i] = 10;
	}

	if (particle) {
		TVector3 True_mom(particle->get_px(), particle->get_py(),
				particle->get_pz());

		seed_mom.SetXYZ(particle->get_px(), particle->get_py(),
				particle->get_pz());
		if(do_smearing)
		{
			const double momSmear = 3. / 180. * TMath::Pi();     // rad
			const double momMagSmear = 0.1;   // relative

			seed_mom.SetPhi(gRandom->Gaus(True_mom.Phi(), momSmear));
			seed_mom.SetTheta(gRandom->Gaus(True_mom.Theta(), momSmear));
			seed_mom.SetMag(
					gRandom->Gaus(True_mom.Mag(),
							momMagSmear * True_mom.Mag()));
		}
	}
}

int get_ioctant(double x, double y) {
	int ioctant = WILD_NUM
	;

	if (x >= 0 && y >= 0) {
		if (y > sqrt(3.) * x)
			ioctant = 0;
		else if (sqrt(3.) * y > x)
			ioctant = 1;
		else
			ioctant = 2;
	} else if (x >= 0 && y < 0) {
		if (sqrt(3.) * -y < x)
			ioctant = 2;
		else if (-y < sqrt(3.) * x)
			ioctant = 3;
		else
			ioctant = 4;
	} else if (x < 0 && y < 0) {
		if (-y > sqrt(3.) * -x)
			ioctant = 4;
		else if (sqrt(3.) * -y > -x)
			ioctant = 5;
		else
			ioctant = 6;
	} else if (x < 0 && y >= 0) {
		if (sqrt(3.) * y < -x)
			ioctant = 6;
		else if (y < sqrt(3.) * -x)
			ioctant = 7;
		else
			ioctant = 0;
	}

	return ioctant;
}

double etaToAngle(double eta) {
	return 2 * TMath::ATan(TMath::Exp(-eta));
}


double quadratic_sum(double a, double b)
{
	return sqrt(a*a + b*b);
}

PHGenFit::PlanarMeasurement* VertexMeasurement(TVector3 vtx, double dr,
		double dphi) {
	PHGenFit::PlanarMeasurement* meas = NULL;

	TVector3 u(1, 0, 0);
	TVector3 v(0, 1, 0);

	TVector3 pos = vtx;

	double u_smear = gRandom->Gaus(0, dphi);
	double v_smear = gRandom->Gaus(0, dr);
	pos.SetX(vtx.X() + u_smear * u.X() + v_smear * v.X());
	pos.SetY(vtx.Y() + u_smear * u.Y() + v_smear * v.Y());

	meas = new PHGenFit::PlanarMeasurement(pos, u, v, dr, dphi);

	return meas;
}

PHGenFit::PlanarMeasurement* PHG4HitToMeasurementVerticalPlane(PHG4Hit* g4hit, int istation =
		0) {
	PHGenFit::PlanarMeasurement* meas = NULL;

	TVector3 pos(
			g4hit->get_avg_x(),
			g4hit->get_avg_y(),
			g4hit->get_avg_z());


	TVector3 v(pos.X(), pos.Y(), 0);
	v = 1/v.Mag() * v;

	TVector3 u = v.Cross(TVector3(0,0,1));
	u = 1/u.Mag() * u;

	double u_smear = gRandom->Gaus(0, _FGEM_phi_resolution);
	double v_smear = gRandom->Gaus(0, _FGEM_r_resolution);
	pos.SetX(g4hit->get_avg_x() + u_smear*u.X() + v_smear*v.X());
	pos.SetY(g4hit->get_avg_y() + u_smear*u.Y() + v_smear*v.Y());


	meas = new PHGenFit::PlanarMeasurement(pos, u, v, _FGEM_phi_resolution, _FGEM_r_resolution);

//	std::cout<<"------------\n";
//	pos.Print();
//	std::cout<<"at "<<istation<<" station, "<<ioctant << " octant \n";
//	u.Print();
//	v.Print();

	//dynamic_cast<PHGenFit::PlanarMeasurement*> (meas)->getMeasurement()->Print();

	return meas;
}


PHGenFit::PlanarMeasurement* PHG4HitToMeasurement(PHG4Hit* g4hit, int istation =
		0) {
	PHGenFit::PlanarMeasurement* meas = NULL;

	TVector3 u(1, 0, 0);
	TVector3 v(0, 1, 0);

	TVector3 pos(
			g4hit->get_avg_x(),
			g4hit->get_avg_y(),
			g4hit->get_avg_z());

	int ioctant = get_ioctant(pos.X(), pos.Y());
	u.RotateZ(-TMath::Pi() / 4. * ioctant);

	double eta = pos.Eta();

	if (istation == 2) {
		v.RotateX(-0.1);
	} else if (istation == 3 || istation == 4) {
		if (eta < 2)
			v.RotateX(-0.1);
		else
			v.RotateX(-0.5 * (etaToAngle(1.45) + etaToAngle(2)));
	}

	v.RotateZ(-TMath::Pi() / 4. * ioctant);

	double u_smear = gRandom->Gaus(0, _FGEM_phi_resolution);
	double v_smear = gRandom->Gaus(0, _FGEM_r_resolution);

	pos.SetX(g4hit->get_avg_x() + u_smear*u.X() + v_smear*v.X());
	pos.SetY(g4hit->get_avg_y() + u_smear*u.Y() + v_smear*v.Y());

	meas = new PHGenFit::PlanarMeasurement(pos, u, v, _FGEM_phi_resolution, _FGEM_r_resolution);

//	std::cout<<"------------\n";
//	pos.Print();
//	std::cout<<"at "<<istation<<" station, "<<ioctant << " octant \n";
//	u.Print();
//	v.Print();

	//dynamic_cast<PHGenFit::PlanarMeasurement*> (meas)->getMeasurement()->Print();

	return meas;
}

double calc_dist_P_to_OA(const TVector3 &P, const TVector3 &OA)
{
	// Project OP to direction that perpenticular to OA -> PC
	TVector3 PC = 1./OA.Mag() * P.Cross(OA);

	// Project PC to direction that perpenticular to OP -> CD; to get PC in phi direction projection
	//TVector3 CD = 1./P.Mag() * PC.Cross(P);

	return PC.Mag();
}

double calc_dist_P_to_AB(const TVector3 &P, const TVector3 &A, const TVector3 &B)
{
	TVector3 AB = B - A;
	return calc_dist_P_to_OA(P, AB) - calc_dist_P_to_OA(A, AB);
}

double calc_dist_P_to_OA_phi(const TVector3 &P, const TVector3 &OA)
{
	// Project OP to direction that perpenticular to OA -> PC
	TVector3 PC = 1./OA.Mag() * P.Cross(OA);

	// Project PC to direction that perpenticular to OP -> CD; to get PC in phi direction projection
	TVector3 CD = 1./P.Mag() * PC.Cross(P);

	return CD.Mag();
}

double calc_dist_P_to_AB_phi(const TVector3 &P, const TVector3 &A, const TVector3 &B)
{
	TVector3 AB = B - A;
	return calc_dist_P_to_OA_phi(P, AB) - calc_dist_P_to_OA_phi(A, AB);
}

int main(int argc, char**argv) {

	if(argc>2)
	{
		_FGEM_phi_resolution = atof(argv[2]);
		LogPrint(_FGEM_phi_resolution);
	}

	if(argc>3)
	{
		_FGEM_r_resolution = atof(argv[3]);
		LogPrint(_FGEM_r_resolution);
	}

	TFile *fout = TFile::Open("out.root", "recreate");
	TTree *Tout = new TTree("T", "Eval NTuple");

	/*!
	 *  Define NTuple vars
	 */

	double Truth_p;
	double Truth_pT;
	double Truth_eta;

	double GenFit_p;
	double GenFit_pT;

	double DCAr;
	double sagitta_st2;

	Tout->Branch("Truth_p", &Truth_p, "Truth_p/D");
	Tout->Branch("Truth_pT", &Truth_pT, "Truth_pT/D");
	Tout->Branch("Truth_eta", &Truth_eta, "Truth_eta/D");


	Tout->Branch("GenFit_p", &GenFit_p, "GenFit_p/D");
	Tout->Branch("GenFit_pT", &GenFit_pT, "GenFit_pT/D");

	Tout->Branch("DCAr", &DCAr, "DCAr/D");
	Tout->Branch("sagitta_st2", &sagitta_st2, "sagitta_st2/D");


	TFile *fPHG4Hits = TFile::Open("G4fsPHENIX.root_DSTReader.root", "read");
	if (!fPHG4Hits) {
		std::cout << "No TFile Openned: " << __LINE__ << "\n";
		return -1;
	}
	TTree *T = (TTree*) fPHG4Hits->Get("T");
	if (!T) {
		std::cout << "No TTree Found: " << __LINE__ << "\n";
		return -1;
	}

	//double nentries = 10000;
	double nentries = T->GetEntries();
	if(argc>1 && atof(argv[1])<nentries)
		nentries = atof(argv[1]);

	//! Initiallize Geometry, Field, Fitter
	bool do_eve_display = false;
	if(nentries<=10) do_eve_display = true;
	const double magnet_field_rescale = 1.0; // 1.0 - 0.0037
	PHGenFit::Fitter* fitter = new PHGenFit::Fitter("fsphenix.root",
			"field.2d.root", magnet_field_rescale, "KalmanFitterRefTrack", "RKTrackRep",
			do_eve_display);


	TH2D *hp_residual_vs_eta = new TH2D("hp_residual_vs_eta",
			"#Delta p/p/p; eta; #Delta p/p/p", 26, 1.45, 4.05, 1000, -0.1, 0.1);

	TH2D *hp_residual_vs_p = new TH2D("hp_residual_vs_p",
			"#Delta p/p; p[GeV/c]; #Delta p/p", 30, 0, 30, 1000, -1, 1);

	TH2D *hsagitta_vs_p = new TH2D("hsagitta_vs_p",
			"#frac{#Delta Sagitta}{Sagitta} (St2) vs. momentum; p[GeV/c]; #frac{#Delta Sagitta}{Sagitta} (St2)", 30, 0, 30, 1000, -1, 1);

	TH2D *hDCAr_vs_p = new TH2D("hDCAr_vs_p",
			"DCAr vs. p; p [GeV/c]; DCAr [cm]", 30, 0, 30, 1000, -1, 1);

#define NLAYERS 5

	Int_t n_G4HIT_FGEM[NLAYERS];
	TClonesArray *G4HIT_FGEM[NLAYERS];
	for (int i = 0; i < NLAYERS; i++) {
		G4HIT_FGEM[i] = NULL;
	}

	T->SetBranchAddress("n_G4HIT_FGEM_0", &n_G4HIT_FGEM[0]);
	T->SetBranchAddress("n_G4HIT_FGEM_1", &n_G4HIT_FGEM[1]);
	T->SetBranchAddress("n_G4HIT_FGEM_2", &n_G4HIT_FGEM[2]);
	T->SetBranchAddress("n_G4HIT_FGEM_3", &n_G4HIT_FGEM[3]);
	T->SetBranchAddress("n_G4HIT_FGEM_4", &n_G4HIT_FGEM[4]);

	T->SetBranchAddress("G4HIT_FGEM_0", &G4HIT_FGEM[0]);
	T->SetBranchAddress("G4HIT_FGEM_1", &G4HIT_FGEM[1]);
	T->SetBranchAddress("G4HIT_FGEM_2", &G4HIT_FGEM[2]);
	T->SetBranchAddress("G4HIT_FGEM_3", &G4HIT_FGEM[3]);
	T->SetBranchAddress("G4HIT_FGEM_4", &G4HIT_FGEM[4]);

	Int_t n_PHG4Particle;
	TClonesArray* TruthParticle = NULL;

	T->SetBranchAddress("n_PHG4Particle", &n_PHG4Particle);
	T->SetBranchAddress("PHG4Particle", &TruthParticle);


	int nevents_nhit_less_than_3 = 0;
	int nevents_failed_fitter = 0;


	for (unsigned int ientry = 0; ientry < nentries; ++ientry) {
		//T->GetEntry(atoi(argv[1]));
		if (ientry % 100 == 0)
			std::cout << "Processing: " << 100. * ientry / nentries << "%"
					<< "\n";

		T->GetEntry(ientry);

		/*!
		 *  Initialize NTuple vars
		 */
		GenFit_p = WILD_NUM;
		Truth_p = WILD_NUM;
		GenFit_pT = WILD_NUM;
		Truth_pT = WILD_NUM;
		DCAr = WILD_NUM;
		sagitta_st2 = WILD_NUM;

		TVector3 seed_pos(0, 0, 0);
		TVector3 seed_mom(0, 0, 0);
		TMatrixDSym seed_cov(6);

		PHG4Particle* particle = NULL;
		for (int i = 0; i < n_PHG4Particle; i++) {
			//LogDEBUG;
			PHG4Particle* temp = dynamic_cast<PHG4Particle*>(TruthParticle->At(
					i));
			//temp->identify();
			//LogDEBUG;
			if (temp->get_track_id() == 1) {
				particle = temp;
				break;
			}
		}

		if (!particle) {
			//LogDEBUG;
			continue;
		}

		get_seed(seed_pos, seed_mom, seed_cov, particle);

		//! Build TrackRep from particle assumption
		/*!
		 * mu+:	-13
		 * mu-:	13
		 * pi+:	211
		 * pi-:	-211
		 * e-:	11
		 * e+:	-11
		 */
		int pid = 13; //
		//SMART(genfit::AbsTrackRep) rep = NEW(genfit::RKTrackRep)(pid);
		genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pid);

		//! Initiallize track with seed from pattern recognition
		PHGenFit::Track* track = new PHGenFit::Track(rep, seed_pos, seed_mom,
				seed_cov);

		TVector3 Hit_St[NLAYERS];
		for(int ilayer = 0; ilayer < NLAYERS; ilayer++) {
			Hit_St[ilayer].SetXYZ(0,0,0);
		}

		//! Create measurements
		std::vector<PHGenFit::Measurement*> measurements;

		const bool use_vertex_in_fitting = true;

		PHGenFit::Measurement* vtx_meas = NULL;

		if (use_vertex_in_fitting) {
			vtx_meas = VertexMeasurement(
					TVector3(0, 0, 0), 0.0050, 0.0050);
			measurements.push_back(vtx_meas);
		}

		for (int ilayer = 0; ilayer < NLAYERS; ilayer++) {
			PHG4Hit* hit = NULL;
			for (int ihit = 0; ihit < n_G4HIT_FGEM[ilayer]; ihit++) {
				PHG4Hit* temp = dynamic_cast<PHG4Hit*>(G4HIT_FGEM[ilayer]->At(
						ihit));

				/*!
				 * TODO for single track simulation, trkid = 1 is the one.
				 * So this code ONLY works for single track simulation right now.
				 */
				if (temp->get_trkid() == 1) {
					hit = temp;
					//hit->identify();
					break;
				}
			}

			if (hit) {
				//PHGenFit::Measurement* meas = PHG4HitToMeasurement(hit, ilayer);
				PHGenFit::Measurement* meas = PHG4HitToMeasurementVerticalPlane(hit, ilayer);

				//meas->getMeasurement()->Print();
				measurements.push_back(meas);

				Hit_St[ilayer].SetXYZ(hit->get_avg_x(),hit->get_avg_y(),hit->get_avg_z());

//				std::cout<<"---------------\n";
//				LogPrint(ilayer);
//				LogPrint(Hit_St[ilayer].Mag());
			}


		}

		//std::cout<<"measurements.size() "<<measurements.size()<<" \n";

		if (measurements.size() < 3) {
			//LogDEBUG;
			nevents_nhit_less_than_3++;
			continue;
		}

		//LogDEBUG;
		//! Add measurements to track
		track->addMeasurements(measurements);

		//LogDEBUG;
		//! Fit the track
		int fitting_err = fitter->processTrack(track, do_eve_display);

		if(fitting_err != 0)
		{
//			LogDEBUG;
//			std::cout<<"event: "<<ientry<<"\n";
			nevents_failed_fitter++;
			continue;
		}
		//LogDEBUG;
		//!
		genfit::MeasuredStateOnPlane* state_at_beam_line =
				track->extrapolateToPlane(TVector3(0, 0, 0), TVector3(0, 0, 1));
		//state_at_beam_line->Print();

		if(!state_at_beam_line)
		{
			LogDEBUG("");
			std::cout<<"event: "<<ientry<<"\n";
			continue;
		}

		TVector3 GenFit_mom = state_at_beam_line->getMom();

		TVector3 True_mom(particle->get_px(), particle->get_py(),
				particle->get_pz());


		Truth_p = True_mom.Mag();
		Truth_pT = True_mom.Pt();
		Truth_eta = True_mom.Eta();

		GenFit_p = GenFit_mom.Mag();
		GenFit_pT = GenFit_mom.Pt();

		//LogDEBUG;
		//hp_residual_vs_p->Fill(True_mom.Pt(),(GenFit_mom.Pt() - True_mom.Pt()) / True_mom.Pt());
		hp_residual_vs_p->Fill(Truth_p, (GenFit_p - Truth_p) / Truth_p);



		hp_residual_vs_eta->Fill(Truth_eta, (GenFit_p - Truth_p) / Truth_p / Truth_p);

		double pu = state_at_beam_line->getState()[1];
		double pv = state_at_beam_line->getState()[2];
		double u = state_at_beam_line->getState()[3];
		double v = state_at_beam_line->getState()[4];
		DCAr = (u*pu + v*pv) / sqrt(pu*pu + pv*pv);

		//LogDEBUG;
		//hDCAr_vs_pT->Fill(True_mom.Mag(), state_at_beam_line->getState()[3]);
		hDCAr_vs_p->Fill(Truth_p, DCAr);

		/*!
		 *  Sagitta at Station 2 calculation, based on GEANT simulation
		 */
//		LogPrint(Hit_St[0].Mag());
//		LogPrint(Hit_St[2].Mag());
//		LogPrint(Hit_St[4].Mag());
		if(
				((use_vertex_in_fitting) || Hit_St[0].Mag() > 0 || Hit_St[1].Mag() > 0) &&
				Hit_St[2].Mag() > 0 &&
				Hit_St[4].Mag() > 0)
		{
			if(use_vertex_in_fitting)
			{
				sagitta_st2 = calc_dist_P_to_AB(Hit_St[2], TVector3(0, 0, 0), Hit_St[4]);
			}
			else if(Hit_St[0].Mag() > 0)
				sagitta_st2 = calc_dist_P_to_AB(Hit_St[2], Hit_St[0], Hit_St[4]);
			else
				sagitta_st2 = calc_dist_P_to_AB(Hit_St[2], Hit_St[1], Hit_St[4]);

			hsagitta_vs_p->Fill(Truth_p, sqrt(1.5)*_FGEM_phi_resolution / sagitta_st2);
		}

		Tout->Fill();

		delete state_at_beam_line;
		delete track;
		//measurements.clear();
		//LogDEBUG;
	} // Loop of events

	std::cout<<"---------------\n";
	LogPrint(nevents_nhit_less_than_3);
	LogPrint(nevents_failed_fitter);
	std::cout<<"---------------\n";

	gStyle->SetOptFit();
	gStyle->SetOptStat(000000000);


	TCanvas *c2 = new TCanvas("c2", "c2");


	hp_residual_vs_eta->FitSlicesY();
	TH1D *hp_resolution_vs_eta = (TH1D*) gDirectory->Get("hp_residual_vs_eta_2");
	hp_resolution_vs_eta->SetTitle(
			"PHGenFit: #sigma_{p}/p/p; #eta; #sigma_{p}/p/p");
	hp_resolution_vs_eta->SetMarkerStyle(20);
	hp_resolution_vs_eta->SetMarkerColor(kBlack);
	hp_resolution_vs_eta->SetLineColor(kBlack);
	hp_resolution_vs_eta->Draw("e");
	hp_resolution_vs_eta->GetYaxis()->SetTitleOffset(1.5);


	TCanvas *c3 = new TCanvas("c3", "c3");
	c3->Divide(2,1);

	c3->cd(1)->SetGrid();
	hp_residual_vs_p->FitSlicesY();
	TH1D *hp_resolution_vs_p = (TH1D*) gDirectory->Get("hp_residual_vs_p_2");
	hp_resolution_vs_p->SetTitle(
			"PHGenFit: #sigma_{p}/p; p[GeV/c]; #sigma_{p}/p");
	hp_resolution_vs_p->SetMarkerStyle(20);
	hp_resolution_vs_p->SetMarkerColor(kBlack);
	hp_resolution_vs_p->SetLineColor(kBlack);
	hp_resolution_vs_p->Draw("e");
	hp_resolution_vs_p->GetYaxis()->SetTitleOffset(1.5);

	TF1 *tf_pT_resolution = new TF1("tf_pT_resolution",
			"sqrt([0]*[0] + x*x*[1]*[1])", 0, 40);
	tf_pT_resolution->SetParameters(0, 0);
	TFitResultPtr fit_result = hp_resolution_vs_p->Fit(tf_pT_resolution,"s");
//	double genfit_mom_resolution_chi2 = fit_result.Get()->Chi2();
//	int genfit_mom_resolution_ndf = fit_result.Get()->Ndf();
	double genfit_mom_resolution_p0 = fit_result.Get()->GetParams()[0];
	double genfit_mom_resolution_p1 = fit_result.Get()->GetParams()[1];

	TH1D *hp_offset_vs_p = (TH1D*) gDirectory->Get("hp_residual_vs_p_1");
	//hp_resolution_vs_p->SetMinimum(-0.01);
	hp_offset_vs_p->GetYaxis()->SetRangeUser(-0.02,0.03);
	hp_offset_vs_p->Draw("e,same");
	hp_offset_vs_p->SetMarkerStyle(20);
	hp_offset_vs_p->SetMarkerColor(kGray);
	hp_offset_vs_p->SetLineColor(kGray);
	fit_result = hp_offset_vs_p->Fit("pol0","s");
	double genfit_mom_offset_p0 = fit_result.Get()->GetParams()[0];

	TH1D *hsagitta_resolution_vs_p = hsagitta_vs_p->ProfileX()->ProjectionX("hsagitta_resolution_vs_p");
	hsagitta_resolution_vs_p->Draw("same");
	hsagitta_resolution_vs_p->SetMarkerStyle(4);
	hsagitta_resolution_vs_p->SetMarkerColor(kBlue);
	hsagitta_resolution_vs_p->SetLineColor(kBlue);
	fit_result = hsagitta_resolution_vs_p->Fit("pol1","s");
	double sagitta_mom_resolution_p0 = fit_result.Get()->GetParams()[0];
	double sagitta_mom_resolution_p1 = fit_result.Get()->GetParams()[1];


	TPaveText *pt = new TPaveText(0.1,0.6,0.7,0.9,"blNDC");
	pt->SetBorderSize(1);
	pt->SetFillColor(0);
	pt->AddText(Form("#frac{#sigma_{p}}{p} = #sqrt{ %f^{2} + (%f*p)^{2}}",genfit_mom_resolution_p0,genfit_mom_resolution_p1));
	pt->AddText(Form("Momentum offset: %f", genfit_mom_offset_p0));
	pt->AddText(Form("#frac{#sigma_{Sagitta}}{Sagitta} = %f + %f*p", sagitta_mom_resolution_p0, sagitta_mom_resolution_p1));

	pt->Draw();

	c3->cd(2);
	hp_residual_vs_p->Draw("colz");
	TH1D *hp_residual_vs_p_mean_sigma = (TH1D*) hp_offset_vs_p->Clone("hp_residual_vs_p_mean_sigma");
	for(int ibin = 1; ibin < hp_residual_vs_p_mean_sigma->GetXaxis()->GetNbins(); ibin++)
	{
		hp_residual_vs_p_mean_sigma->SetBinError(ibin, hp_resolution_vs_p->GetBinContent(ibin));
	}
	hp_residual_vs_p_mean_sigma->Draw("same");
	hp_residual_vs_p_mean_sigma->SetMarkerStyle(4);
	hp_residual_vs_p_mean_sigma->SetMarkerColor(kRed);
	hp_residual_vs_p_mean_sigma->SetLineColor(kRed);



	TCanvas *c4 = new TCanvas("c4", "c4");


	hDCAr_vs_p->FitSlicesY();
	TH1D *hDCAr_resolution_vs_p = (TH1D*) gDirectory->Get("hDCAr_vs_p_2");
	hDCAr_resolution_vs_p->SetTitle(
			"PHGenFit: #sigma_{DCAr}; p [GeV/c]; #sigma_{DCAr} [cm]");
	hDCAr_resolution_vs_p->SetMarkerStyle(20);
	hDCAr_resolution_vs_p->SetMarkerColor(kBlack);
	hDCAr_resolution_vs_p->SetLineColor(kBlack);
	hDCAr_resolution_vs_p->Draw("e");
	hDCAr_resolution_vs_p->GetYaxis()->SetTitleOffset(1.5);
	//hDCAr_resolution_vs_p->SetMinimum(-0.01);
	hDCAr_resolution_vs_p->GetYaxis()->SetRangeUser(-0.02,0.1);


	TF1 *tf_DCAr_resolution = new TF1("tf_DCAr_resolution",
				"sqrt([0]*[0] + [1]*[1]/x/x)", 0, 40);
	tf_DCAr_resolution->SetParameters(0,0);
	hDCAr_resolution_vs_p->Fit(tf_DCAr_resolution);

	TH1D *hDCAr_offset_vs_p = (TH1D*) gDirectory->Get("hDCAr_vs_p_1");

	hDCAr_offset_vs_p->Draw("e,same");
	hDCAr_offset_vs_p->SetMarkerStyle(20);
	hDCAr_offset_vs_p->SetMarkerColor(kGray);
	hDCAr_offset_vs_p->SetLineColor(kGray);
	hDCAr_offset_vs_p->Fit("pol0");

	TCanvas *c0 = new TCanvas("c0", "c0");
	hsagitta_vs_p->Draw();
	hsagitta_resolution_vs_p->Draw("same");
	c0->Update();

	//! Event display
	if (do_eve_display)
	{
		fitter->displayEvent();
		std::cin.get();
	}

	fPHG4Hits->Close();

	fout->cd();
	Tout->Write();
	c0->Write();
	c2->Write();
	c3->Write();
	hp_resolution_vs_p->Write();
	c4->Write();
	fout->Close();

	//delete fitter;

	//pause();

	std::cout << "SUCCESS! \n";

	return 0;
}
